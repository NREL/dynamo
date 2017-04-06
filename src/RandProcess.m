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
%  12  2017-04-05 23:52  BryanP      Added checkState, varios bugfixes
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

    % Read only properties
    properties (GetAccess = 'public', SetAccess='protected')
        t = NaN           %current timestep
        cur_state = NaN   %current state
    end

    methods (Abstract)
        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.
        %
        % IMPORTANT: all processes must have a single, starting state for t=0
        % DLIST List possible discrete states
        %
        % List possible discrete states by number for given time
        % if t is not listed, the states for the current simulation time
        % are returned.
        %
        % To get a list of all possible states pass with t='all', in which
        % case, t_max is provided to bound the problem size.
        state_list = dlist (obj, t)

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

        % DSIM Simulate discrete process.
        %
        % A column vector for t is assumed to be a series of times for
        % which to return results. Intermediate times are also computed, if
        % needed, but not returned. The initial value is not returned in
        % the value series. Only one simulation is run, such that out of
        % order times will be sorted before simulation and duplicate times
        % will return the same result
        %
        % Invalid times (t<1) return NaN
        %
        % Note: after calling sim, the process internal time will be set to
        % the final value of t_list
        state_series = dsim(obj, t_list, initial_value)

        %% ===== General (discrete or continuous) Methods

        % SIM Simulate process for desired (continuous) times
        %
        % A column vector for t is assumed to be a series of times for
        % which to return results. Intermediate times are also computed, if
        % needed. The initial value is not returned in the value series.
        %
        % Function must handle arbitrary positive values for t_list
        % Invalid times return NaN
        %
        % Note: after calling sim, the process internal time will be set to
        % the final value of t_list
        state_series = sim(obj, t_list, initial_value)

        % RANGE Find value range for given time
        %
        % Returns vector with [min max] value range for specified time
        % if t is not provided, the range for the current simulation time
        % is returned.
        %
        % To get the possible range across all times use t='all'
        state_value_range = range(obj, t)
        
        %SAMPLE draw state samples for the given time
        %
        % Usage:
        %   state_list = disc_samp_object.sample()
        %       One sample state from current time
        %   state_list = sample(obj, N)
        %       Return N samples from current time
        %   state_list = sample(obj, N, t)
        %       Specify time period
        state_list = sample(obj, N, t)

        %% ===== Additional simulation support
        %STEP simulate forward
        %
        % by default steps forward by delta_t = 1
        [state, t] = step(obj, delta_t)
        
        % CHECKSTATE Check that state is valid for a given time
        %
        % rand_proc_objectcCheckState(t, state)
        %       Raise 'RandProc:InvalidState' error if t is not valid in time t
        % state_ok rand_proc_object.checkState(t, state)
        %       No error, simply return true if state is valid
        state_ok = checkState(obj, t, state)

    end

    methods

        function [val, prob] = as_array(obj, varargin)
            %For a RandProcess, as_array is a simple wrapper around
            %dlistnext for the current state and time.c

            [val, prob] = obj.dlistnext(obj.cur_state, obj.t, varargin);
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
        
    end

end
