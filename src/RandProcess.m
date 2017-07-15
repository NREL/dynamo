classdef RandProcess < AbstractSet
%RANDPROCESS Random Process abstract class for dynamic programming
%
% Defines an abstract class (defines the strucuture for related
% subclasses) for a Random Process for use with dynamic programming and
% other simulations
%
% IMPORTANT: all processes must have a single, starting state for t=1
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  16  2017-07-14 21:35  BryanP      Remove sim and dsim, distinquish N_uniqueT from Tmax 
%  15  2017-07-14 11:10  BryanP      Refactor generally usable code from rpDiscreteSample to here 
%  14  2017-07-14 06:02  BryanP      Clarify condintional and unconditional use for sample() 
%  13  2017-04-10 16:36  BryanP      Make sample depend on state 
%  12  2017-04-05 23:52  BryanP      Added checkState, various bugfixes
%  11  2017-04-05 22:52  BryanP      include reset() implementation
%  10  2017-04-05 22:02  BryanP      Put t as first parameter for dlistnext
%   9  2017-04-05 15:12  BryanP      Setup for 1-based t indexing and pure value states (no state_n)  
%   8  2016-10-07 01:40  BryanP      Convert as_array to wrapper for dlistnext 
%   7  2016-10-06 11:40  DheepakK    as_array that throws an error
%   6  2016-07-07 09:40  BryanP      Extracted out common pieces to form AbstractSet class
%   5  2010-12-23 19:30  BryanP      Made t an abstract property to allow subclass error handling
%   4  2010-12-15 23:30  BryanP      Added t_max for dlist(... 'all')
%   3  2010-12-14 12:44  BryanP      Added value returns to dlist*
%   2  2010-12-13 21:20  BryanP      distinguish dsim & sim add internal t
%   1  2010-12-13 12:20  BryanP      Initial Version

    % User modifiable tolerances
    properties
        Tol = 1e-6;     % Tolerance for checking probabilities and state membership
    end
    
    % Read only properties
    properties (GetAccess = 'public', SetAccess='protected')
        t = NaN           %current timestep
        cur_state = NaN   %current state
    end
    
    % Internal properties
    properties (Access='protected')
        Values = {};    % cell row vector of values: 1 cell column per time, typically column vectors per time
        UncondProbs = {};      % cell row vector of (unconditional) probabilities for each time period
        UncondCdfs = {}; % cell row vector of (unconditional) cumulative distribution for each time period
        
        Tmax = Inf;       % Largest time with specified distribution
    end

    properties (Dependent, Access = protected)
        N_uniqueT;   % Number of unique time periods of values defined
    end

    methods (Abstract)
        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.
        %
        % IMPORTANT: all processes must have a single, starting state for
        % t=1

        % DLISTNEXT List next discrete states & probabilities
        %
        % List possible next states (by number) along with conditional
        % probability P(s_{t+1} | s_t)
        %
        % If t is not provided, the current simulation time is assumed
        %
        % If the state is not valid at t, an error with ID
        % 'RandProcess:InvalidState' should be thrown
        %
        % If t is out of the valid range, an error with ID
        % 'RandProcess:InvalidTime'
        [next_state_list, prob] = dlistnext (obj, t, state )
    end
    
    methods (Abstract, Access = protected)
        %% ===== Process specific sub-functions
        %CONDITIONALSAMPLE draw state samples for the specified state
        %
        % This helper function implements the conditional probability
        % support specific to the subclass
        %
        % CONDITIONAL PROBABLITITY SUPPORT (Implemented in subclass via conditionalSample method)
        %   state_list = sample(obj, N, t, cur_state)
        %       Sample specified time using conditional probability
        %       starting from provided state
        %
        %   Set cur_state to empty to use the object's current state
        state_list = conditionalSample(obj, N, t, cur_state)
    end

    %% =========== General Public methods
    methods
        function state_list = sample(obj, N, t, cur_state)
        %SAMPLE draw state samples for the given time and (current) state
        %
        % Usage:
        %   state_list = disc_samp_object.sample()
        %       One sample state from current time, using conditional
        %       probability for current_state
        %   state_list = sample(obj, N)
        %       Return N samples
        %   state_list = sample(obj, N, t)
        %       Specify time and sample based on unconditional probability
        %       across all valid states for t
        %
        % CONDITIONAL PROBABLITITY SUPPORT (Implemented in subclass via conditionalSample method)
        %   state_list = sample(obj, N, t, cur_state)
        %       Sample specified time using conditional probability
        %       starting from provided state
            if nargin < 2
                N = 1;
            end
            
            if nargin < 3 || isempty(t)
                t=obj.t;
            end
            
            if nargin >= 4 && not(isempty(cur_state))
                state_list = obj.conditionalSample(N, t, cur_state);
                return
            end
            
            %Handle any non integer or large values for time
            t = min(floor(t), obj.N_uniqueT);

            idx_at_t = zeros(N,1);
            for samp_idx = 1:N
                idx_at_t(samp_idx) = find(rand(1) <= obj.UncondCdfs{t}, 1, 'first');
            end
            state_list = obj.Values{t}(idx_at_t,:);
        end


        function [val, prob] = as_array(obj, t, varargin)
        %AS_ARRAY return possible values and (unconditional) probabilities
        % val = rand_proc_object.as_array()
        %   Return list of values for the current time
        %
        % [val, prob] = rand_proc_object.as_array()
        %   Also return the corresponding (unconditional) probability
        %
        % __ = rand_proc_object.as_array(t)
        %   Specify the time to use for the list. (Does not update the
        %   object's current time)
        %
        % val = rand_proc_object.as_array('all')
        %   List all possible states for all times (probability undefined

            
            if nargin < 2 || isempty (t)
                t = obj.t;
            end

            if (ischar(t) && strcmp(t, 'all'))
                val = unique(cell2mat(obj.Values'),'rows');
                prob = NaN;
                return
            elseif t > obj.Tmax
                t = obj.Tmax;
            end

            val = obj.Values{t};
            prob = obj.UncondProbs{t};
        end

        function reset(obj, initial_state)
            % RESET reset simulation
            %
            % rand_proc_obj.reset()
            %       resets t=1 and a random initial state
            % rand_proc_obj.reset(initial_state)
        
            obj.t = 1;
            if nargin > 1
                %Check state will error out if state is invalid
                obj.checkState(obj.t, initial_state);
                
                obj.cur_state = initial_state;
            else
                obj.cur_state = obj.sample();
            end
        end

        function value_range = range(obj, t)
        % RANGE Find value range for given time
        %
        % Returns vector with [min max] value range for specified time
        % if t is not provided, the range for the current simulation time
        % is returned.
        %
        % To get the possible range across all times use t='all'
            if nargin < 2 || isempty(t)
                t= obj.t;
            end

            if ischar(t) && strcmp(t, 'all')
                state_list_to_range = cell2mat(obj.Values');
            else
                %Handle any non-integer or large values
                t = min(floor(t), obj.N_uniqueT);

                if t < 1
                    error('RandProcess:InvalidTime', 'Only t>=1 valid for Random Processes')
                else
                    state_list_to_range = obj.Values{t};
                end
            end
            value_range = [min(state_list_to_range); max(state_list_to_range)];
        end


        %% ===== Additional simulation support
        function [state, t] = step(obj, delta_t)
        %STEP simulate forward from current (internal) state
        %
        % by default steps forward by delta_t = 1
            if nargin < 2 || isempty(delta_t)
                delta_t = 1;
            end
            %compute the proposed new time
            new_t = obj.t + delta_t;

            %check if it is valid, if not return empty results
            if new_t < 1 || new_t > obj.Tmax
                state = NaN;
                t = obj.t;
                return
            else
                %if new time is valid simulate forward as needed
                state = obj.cur_state;
                for t = obj.t:new_t-1
                    state = obj.sample(1, t, state);
                end
                % finish by moving forward the last time step if needed
                if delta_t >= 1
                    t = t+1;
                else
                    t = obj.t;
                end
    
                % Update our stored state
                obj.t = t;
                obj.cur_state = state;
            end
        end
        
        function state_ok = checkState(obj, t, state)
        % CHECKSTATE Check that state is valid for a given time
        %
        % rand_proc_object.checkState(t, state)
        %       Raise 'RandProc:InvalidState' error if t is not valid in time t
        % state_ok = rand_proc_object.checkState(t, state)
        %       No error, simply return true/false if state is
        %       valid/not

            %First double check the time has a chance of being valid
            if t < 1
                state_ok = false;
                if nargout == 0
                    error('RandProcess:InvalidTime', 'Random Process times must be >=1 (not %d)', t)
                end
                return
            elseif t > obj.Tmax
                state_ok = false;
                if nargout == 0
                    error('RandProcess:InvalidTime', 'Random Process times must be <Tmax (not %d)', t)
                end
                return
            end
            
            % Convert time into lookup value
            t_lookup = min(floor(t), obj.N_uniqueT);
                
            state_ok = not(isempty(state)) && ismembertol(state, obj.Values{t_lookup}, obj.Tol, 'ByRows', true);
            
            if nargout == 0 && not(state_ok)
                error('RandProcess:InvalidState', 'State [ %s] is not valid at time %d', sprintf('%g ',state), t)
            end
        end

        % === Reference value functions
        function n = get.N_uniqueT(obj)
            n = length(obj.Values);
        end

    end
    
    methods (Access = protected)
        
        %% ===== Helper Functions
        function [state_list, prob] = state_info (obj, t )
            % STATE_INFO Helper function to return full set of state
            % information for a given time (for discrete time processes)
            if t < 1
                error('RandProcess:InvalidTime', 'Only t>1 valid for rpDiscreteSample')
            else
                % make sure time is an integer
                t = floor(t);
                
                % and use t=Tmax for any t>Tmax
                t = min(t, obj.N_uniqueT);
                
                %if we get here, we know the time is valid
                state_list = obj.Values{t};
                prob = obj.UncondProbs{t};
            end
        end
        
    end

end
