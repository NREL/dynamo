classdef rpDeterministic < RandProcess
%rpList Simple random process with deterministic outcome from list
%
% obj = rpList(v_list)
%
% A "random" process that is deterministic, but for which the value is
% determined by a (discrete) list.
%
% The ValueList maybe passed to the constructor or set directly. It is a
% vector of values corresponding to the fixed process outcome for times
% from t=0 to length(ValueList)-1.
%
% Notes:
%  - IMPORTANT: the value list starts at t=0!
%  - For any t greater than the max set in the ValueList, the final value
%    from the list is assumed to be held constant
%  - For non-integer times, the output is assumed to remain at the value
%    corresponding to the previous integer time step (Zero-order hold)
%
% see also RandProc, rpDiscreteSample, rpMarkov, rpList, rpLattice
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-12-13 12:30  BryanP      Adapted from RandProcess v1
%   2  2010-12-13 21:20  BryanP      Renamed rpList, adapted for RandProcess v2
%   3  2010-12-13 23:20  BryanP      Explicit lookup for StateNum to skip find
%   4  2010-12-14 12:44  BryanP      Added value returns to dlist*
%   5  2010-12-14 23:30  BryanP      Corrected based on unit testing with testrpList
%   6  2010-12-14 23:30  BryanP      Converted to first value being t=0
%   7  2010-12-15 23:30  BryanP      Added t_max for dlist(... 'all')
%   8  2010-12-23 19:30  BryanP      Declare t property to match RandProcess v5
%   9  2011-01-04 22:00  BryanP      Errors for out of range t, etc. consistant with rpLattice v5
%  10  2011-06-09 16:00  BryanP      match input shape for dnum2val and dval2num

    properties
        t = 0;        % Current time
        ValueSeries = [];   % List of possible values, with one entry per time period
    end

    % Internal properties
    properties (Access='protected')
        ValueMap = [];       % List of unique values in ValueSeries, use find for lookup
        StateNumSeries = []; % ValueSeries converted to corresponding state numbers
        n = 0;               % Total number of specified values in ValueSeries
    end

    methods
        %% ===== Constructor
        function obj = rpList(v_list)
            % Support zero parameter calls to constructor for special
            % MATLAB situations (see help files)
            if nargin > 0
                obj.ValueSeries = v_list;
            end
        end

        %% ===== Property value maintenance
        % Set associated properties when assigning ValueSeries
        function set.ValueSeries(obj, values)
            obj.ValueSeries = values;
            [obj.ValueMap, junk, obj.StateNumSeries] = unique(obj.ValueSeries); %#ok<ASGLU,MCSUP>
            obj.n = length(obj.ValueSeries); %#ok<MCSUP>
        end

        % Check that we set a valid t
        function set.t(obj, t)
            if t < 0;
                error('RandProcess:InvalidTime', 'time (t) must be positive')
            else
                obj.t = t;
            end
        end


        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.

        function [value_list, state_n_list] = dlist (obj, t)
        % DLIST List possible discrete states
        %
        % List possible discrete states by number for given time
        % if t is not listed, the states for the current simulation time
        % are returned.
        %
        % To get a list of all possible states pass with t='all', in which
        % case, t_max is required to bound the problem size.
            if nargin < 2 || isempty(t)
                [value_list, state_n_list] = obj.dlist(obj.t);
            elseif isempty(obj.ValueSeries) || (ischar(t) && strcmp(t, 'all'))
                value_list = obj.ValueMap;
                state_n_list = 1:length(obj.ValueMap);
            elseif t < 0;
                error('RandProcess:InvalidTime', 'Only t>0 valid for rpList')
            else
                state_n_list = obj.StateNumSeries(min(obj.n,t+1));
                value_list = obj.ValueMap(state_n_list);
            end
        end

        function [value_list, state_n_list, prob] = dlistprev (obj, state_n, t )
        % DLISTPREV List previous discrete states & probabilities
        %
        % List possible previous states (by number) along with conditional
        % probability P(s_t | s_{t-1})
        %
        % If t is not defined, the current simulation time is assumed
            if nargin < 3
                if nargin < 2
                    state_n = obj.StateNumSeries(min(obj.n,obj.t+1));
                end
                [value_list, state_n_list, prob] = obj.dlistprev(state_n, obj.t);
            elseif t <= 0;
                error('RandProcess:InvalidTime', 'Only t>0 valid for rpList')
            elseif isempty(state_n) || state_n == obj.StateNumSeries(min(obj.n,t+1))
                state_n_list = obj.StateNumSeries(min(obj.n,t));
                value_list = obj.ValueMap(state_n_list);
                prob = 1;
            else
                error('RandProcess:InvalidState', 'State #%d not valid at time=%d', state_n, t)
            end
        end

        function [value_list, state_n_list, prob] = dlistnext (obj, state_n, t )
        % DLISTNEXT List next discrete states & probabilities
        %
        % List possible next states (by number) along with conditional
        % probability P(s_{t+1} | s_t)
        %
        % If t and/or state_n are not defined, the current simulation time
        % and state are assumed
            if nargin < 3
                if nargin < 2
                    state_n = obj.StateNumSeries(min(obj.n,obj.t+1));
                end
                [value_list, state_n_list, prob] = obj.dlistnext(state_n, obj.t);
            elseif t < -1;
                error('RandProcess:InvalidTime', 'Only t>0 valid for rpList')
            else
                state_n_list = obj.StateNumSeries(min(obj.n,t+2));
                value_list = obj.ValueMap(state_n_list);
                prob = 1;
            end
        end

        function values = dnum2val (obj, state_n_list )
        % DNUM2VAL Convert discrete state number to a value
            try
                values = obj.ValueMap(state_n_list);
            catch exception
                bad = find(or(state_n_list > obj.n, state_n_list < 1), 1);
                error('rpList:InvalidStateNum', ...
                    'Attempt to find value for invalid state nums: %d', ...
                    state_n_list(bad))
            end
            if size(values) ~= size(state_n_list)
                values = reshape(values, size(state_n_list));
            end
        end

        function state_n_list = dval2num (obj, values )
        % DVAL2NUM Convert the value(s) to associate discrete state numbers
            try
                state_n_list = arrayfun(@(x) find(obj.ValueMap == x, 1), values);
            catch exception
                bad = find(not(ismember(values, obj.ValueMap)), 1);
                error('rpList:InvalidValue', ...
                    'Attempt to find value for invalid state nums: %d', ...
                    values(bad))
            end
            if size(values) ~= size(state_n_list)
                state_n_list = reshape(state_n_list, size(values));
            end
        end

        function [value_series, state_n_series] = dsim(obj, t_list, initial_value) %#ok<INUSD>
        % DSIM Simulate discrete process.
        %
        % A column vector for t is assumed to be a series of times for
        % which to return results. Intermediate times are also computed, if
        % needed, but not returned. The initial value is not returned in
        % the value series. Only one simulation is run, such that out of
        % order times will be sorted before simulation and duplicate times
        % will return the same result
        %
        % Invalid times (t<0) return NaN
        %
        % Note: after calling sim, the process internal time will be set to
        % the final value of t_list
            ok = t_list >= 0;
            value_series = zeros(size(t_list));
            value_series(not(ok)) = NaN;

            %Keep only good values of t_list
            t_list = t_list(ok);
            value_series(ok) = obj.ValueSeries(min(obj.n, t_list+1));
            %if there are valid values from t_list, set t to the last one.
            if not(isempty(t_list))
                obj.t = t_list(end);
            end

            %Handle state numbers if required
            if nargout > 1;
                state_n_series = zeros(size(value_series));
                state_n_series(not(ok)) = NaN;
                state_n_series(ok) = obj.StateNumSeries(min(obj.n, t_list+1));
            end
        end

        %% ===== General (discrete or continuous) Methods
        function [value_series, state_n_series] = sim(obj, t_list, initial_value)
        % SIM Simulate process for desired (continuous) times
        %
        % A column vector for t is assumed to be a series of times for
        % which to return results. Intermediate times are also computed, if
        % needed. The initial value is not returned in the value series.
        %
        % Function must handle arbitrary positive values for t_list
        % Invalid times (t<=0) return NaN
        %
        % Note: after calling sim, the process internal time will be set to
        % the final value of t_list
        %
        % This implementation typically rounds down to the nearest integer
        % (zero order hold)
            if nargin < 3
                initial_value = [];
            end
            if nargout == 1
                value_series = obj.dsim(floor(t_list), initial_value);
            else
                [value_series, state_n_series] = obj.dsim(floor(t_list), initial_value);
            end

        end

        function [value_range, state_n_range] = range(obj, t)
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

            if isempty(obj.ValueSeries) || (ischar(t) && strcmp(t, 'all'))
                value_range = obj.ValueMap([1, end]);
                state_n_range = [1, length(obj.ValueMap)];
            elseif t < 0;
                error('RandProcess:InvalidTime', 'Only t>0 valid for rpList')
            else
                value_range = obj.ValueSeries(min(obj.n,floor(t)+1)) .* [1 1];
                if nargout > 1
                    state_n_range = obj.StateNumSeries(min(obj.n,floor(t)+1)) .* [1 1];
                end
            end
        end


        %% ===== Additional simulation support
        function [value, state_n, t] = step(obj, delta_t)
        %STEP simulate forward
        %
        % by default steps forward by delta_t = 1
            if nargin < 2
                delta_t = 1;
            end
            new_t = obj.t + delta_t;
            obj.t = max(0, new_t);

            if new_t < 0
                value = [];
                state_n = [];
                t = obj.t;
            else
                if nargout == 1
                    value = curState(obj);
                else
                    [value, state_n, t] = curState(obj);
                end
            end
        end

        function [value, state_n, t] = curState(obj)
        %CURSTATE Return the current state of the simulation
            value = obj.ValueSeries(min(obj.n,floor(obj.t+1)));
            if nargout > 1
                state_n = obj.StateNumSeries(min(obj.n,floor(obj.t+1)));
                t = obj.t;
            end
        end


        function reset(obj, initial_value) %#ok<INUSD>
        %RESET reset simulation to t=0
            obj.t = 0;
        end

    end

end
