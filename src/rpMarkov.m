classdef rpMarkov < RandProcess

% WARNING--DO NOT USE. Not updated for latest RandProcess styles
    
%rpMarkov Generalized dynamic Markov random process
%
% Supports a generalized dynamic Markov random process where the
% state values and/or transition probabilities can be a function of time.
%
% Constructor: obj = rpMarkov(v_set, trans_set)
%
% Required inputs/properties:
%  v_set: cell row vector with one cell per specified time period
%           containing a numeric ROW vector of possible values (states).
%  trans_set: matching cell vector with corresponing transition probability
%           matrices between the states at t-1 and t.  Transition matrix
%           entries contain the probability of transitioning from the
%           previous row state to the current column state. Since the
%           length of v_set entries may change between time periods, the
%           matrices may not be square, but rather must have must have as
%           many rows as states in the previous period (ie #rows =
%           length(v_set{t-1})) and as many columns as states in the
%           current period (ie #cols=length(v_set{t})). The first entry
%           provides the initial probability of taking each value from
%           v_set{1}, and therefor should only have one row. The sum of
%           each row must equal 1.0
%
%  It is acceptable for the v_set and trans_set cell arrays to be of
%  different lengths, in which case, the last value from the shorter list
%  will be used. For instance this allows a number of possible Markov
%  derived processes:
%
%   TRADITIONAL MARKOV
%     This corresponds to a traditional Markov chain where the values for
%     each state are constant and the transition probabilies are also
%     constant. To set this up provide the v_set cell array with one entry
%     and trans_set cell array with two elements: the v_set entries should
%     contain a list of n possible values, while the first trans_set entry
%     (trans_set{1}) provides a row vector with the inital probability of
%     being in the corresponding state. For a fixed initial state, specify
%     a one for the corresponding value and zeros for all others. The
%     second trans_set entry then should contain a square n-by-n matrix of
%     transition probabilities.
%   DYNAMIC VALUES
%     This case can be thought of as having a fixed set of abstract states
%     for which the actual value changes over time. For example there might
%     be a low, medium, and high priced scenario, for which the actual
%     prices increase over time. In this case the trans_set should be setup
%     as in a traditional markov, but the v_set should have as many values
%     as required to specify changes over time.
%   DYNAMIC TRANSITIONS
%     In this case the probabilities of transitioning between a fixed set
%     of states changes over time. Setup for v_set is the same as for a
%     traditional markov, but the trans_set cell array should have as many
%     elements are required to capture the changing probabilites.
%   FULLY DYNAMIC
%     Here both values and transitions may change with time.
%
%
% Additional properties (set using obj.PROPERTY):
%  None
%
% Notes:
%  - IMPORTANT: the value list starts at t=0!
%  - For non-integer times, the output is assumed to remain at the value
%    corresponding to the previous integer time step (Zero-order hold)
%
% see also RandProc, rpDiscreteSample, rpMarkov, rpList, rpLattice
%
% originally by Bryan Palmintier 2011

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-05-02 15:30  BryanP      Adapted from rpDiscreteSample v2
%   2  2011-05-03 15:00  BryanP      Initial coding & basic debug
%   3  2011-06-09 16:00  BryanP      match input shape for dnum2val and dval2num
%   4  2012-04-17 08:55  BryanP      Reworked cdf for sim/sample


    properties
        t = 0;        % Current time
    end

    % Read only properties
    properties (GetAccess = 'public', SetAccess='protected')
        VSet = {};    % cell row vector of values
        TransSet = {};    % cell row vector of probabilities

        cur_state = NaN;    %Current state number for simulation
%     end
%
%     % Internal properties
%     properties (Access='protected')
        SSet = {};     % cell row vector of state numbers
        CdfSet = {};   % cumulative distribution for each state in each time period

        ValueMap = [];  % List of unique values in VSet, index provides state number
        IdxMap = cell(0); % Lookup tables to convert state numbers to index numbers for each time period

        nVal = 0;       % number of specifed values
        nTrans = 0;     % number of specifed transition matrices
        Tol = 1e-6;     % Tolerance for checking probabilities sum to 1
    end

    methods (Access = protected)
        %% ===== Helper Functions
        function [state_n, value] = step_helper(obj, s, t)
            %STEP_HELPER step forward one period with no checks

            [v_idx, next_trans] = obj.t2idx(t);

            % Compare a random number on [0,1) to the CDF row corresponding
            % to the current state in the transition matrix
            idx_at_t = find(rand() <= obj.CdfSet{next_trans}(obj.s2idx(s,t),:), 1, 'first');

            % Return is a column number, so need to convert to a state
            % number:
            state_n = obj.SSet{v_idx}(idx_at_t);
            value = obj.VSet{v_idx}(idx_at_t);
        end

        function [idx] = s2idx(obj, s_list, t)
            %S2IDX convert a state numbers to a transition matrix indices
            t = min(floor(t)+1, obj.nVal);

            idx = obj.IdxMap{t}(s_list);
        end

        function [v_idx, cur_trans] = t2idx(obj, t)
            %T2IDX find the value and transition set indexes for time

            if t < 0
                error('RandProcess:InvalidTime', 'Only t>0 valid for rpMarkov')
            else
                % make sure time is an integer
                t = floor(t);

                % and use t=Tmax for any t>Tmax for value and transition
                % matrices independently
                v_idx = min(t + 1, obj.nVal);
                cur_trans = min(t + 1, obj.nTrans);
            end
        end
    end

    methods
        %% ===== Constructor & related
        function obj = rpMarkov(v_set, trans_set)
            % Support zero parameter calls to constructor for special
            % MATLAB situations (see help files)
            if nargin > 0
                obj.VSet = v_set;
                obj.TransSet = trans_set;

                obj.nVal = length(v_set);
                obj.nTrans = length(trans_set);
                obj.ValueMap = unique(vertcat(v_set{:}));

                obj.SSet = cell(size(v_set));
                obj.CdfSet = cell(size(trans_set));
                obj.IdxMap = cell(size(v_set));

                % Loop through time for setting up additional parameters
                % First values
                for t_idx = 1:obj.nVal
                    %Create State list
                    [~, obj.SSet{t_idx}, org_order] = intersect(obj.ValueMap, v_set{t_idx});
                    obj.SSet{t_idx}(org_order) = obj.SSet{t_idx};

                    %And index lookup (result will be a vector with
                    %non-zero entries for all valid states at time
                    %t_idx with values equal to the corresponding indexes
                    %for the value list and transition matrices
                    obj.IdxMap{t_idx} = zeros(size(obj.ValueMap));
                    obj.IdxMap{t_idx}(obj.SSet{t_idx}) = 1:length(obj.SSet{t_idx});
                end

                % Now  handle transition matrices
                for t_idx = 1:obj.nTrans
                    %compute cdf by summing across rows
                    obj.CdfSet{t_idx} = cumsum(obj.TransSet{t_idx}, 2);

                    %Check for valid probability vectors
                    if any(abs(obj.CdfSet{t_idx}(:,end) - 1) > obj.Tol)
                        error('RandProc:ProbSum1', 'Probabilities must sum to 1.0 for t=%d', t_idx-1)
                    end

                    %And valid transition matrix sizes
                    if t_idx == 1
                        %First period should have only one row of same
                        %length as VSet{1}
                        if size(obj.TransSet{1}) ~= size(obj.VSet{1})
                            error('rpMarkov:ValTransMismatch', 'First Period Trans Matrix must be same size as Value row vector')
                        end
                    else
                        %Other periods should have #col = old VSet
                        % and #rows = new/current VSet
                        [r, c] = size(obj.TransSet{t_idx});
                        if r ~= length(obj.VSet{t_idx-1}) || c ~= length(obj.VSet{t_idx})
                            error('rpMarkov:ValTransMismatch', 'Transition Matrix must have #rows=length(v_set(t-1)) and #col=length(v_set(t)) for time=%d', t_idx)
                        end
                    end
                end
                obj.reset();
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
        % To get a list of all possible states pass with t='all'
            if nargin < 2 || isempty(t)
                [value_list, state_n_list] = obj.dlist(obj.t);
            elseif (ischar(t) && strcmp(t, 'all'))
                value_list = obj.ValueMap';
                state_n_list = 1:length(obj.ValueMap);
            else
                v_idx = obj.t2idx(t);
                value_list = obj.VSet{v_idx};
                state_n_list = obj.SSet{v_idx};
            end
        end

        function [value_list, state_n_list, prob] = dlistprev (obj, state_n, t )
        % DLISTPREV List previous discrete states & probabilities
        %
        % List possible previous states (by number) along with conditional
        % probability P(s_{t-1} | s_t). Since the samples are independant
        % of each other (ie not path dependant), this simplifies to
        % P(s_{t-1})
        %
        % If t and/or state_n are not defined, the current simulation time
        % and state are assumed
            if nargin < 3 || isempty(t)
                [value_list, state_n_list, prob] = obj.dlistprev([], obj.t);
                return
            end

            error('RandProcess:NotImplemented', 'The method dlistprev is not yet implemented for rpMarkov')
            %TODO: implement this... will require recursively multiplying
            %the transition probability matrices to reach the current time
            %then applying Bayes theory: doable but a bit of a hassle
        end

        function [value_list, state_n_list, prob] = dlistnext (obj, state_n, t )
        % DLISTNEXT List next discrete states & probabilities
        %
        % List possible next states (by number) along with conditional
        % probability P(s_{t+1} | s_t) Since the samples are independant
        % of each other (ie not path dependant), this simplifies to
        % P(s_{t+1})
        %
        % If t and/or state_n are not defined, the current simulation time
        % and state are assumed
            if nargin < 3
                [value_list, state_n_list, prob] = obj.dlistnext(state_n, obj.t);
                return
            end

            % Get the appropriate set indices for the next time period
            [v_idx, next_trans] = obj.t2idx(t+1);

            % Check for a valid state at this time period
            if nargin > 2 && not(isempty(state_n)) && ...
                        ( state_n > length(obj.ValueMap) ...
                          || not(all(ismember(state_n,obj.SSet{t+1})))...
                        )
                error('RandProcess:InvalidState', 'State #%d is not valid at time %d', state_n, t)
            else
                %if we get here, we know the state is valid so return the
                %required data

                prob = obj.TransSet{next_trans}(obj.s2idx(state_n, t),:);
                value_list = obj.VSet{v_idx}(prob > 0);
                state_n_list = obj.SSet{v_idx}(prob > 0);
                prob = prob(prob > 0);
            end
        end

        function values = dnum2val (obj, state_n_list )
        % DNUM2VAL Convert discrete state number to a value
            try
                values = obj.ValueMap(state_n_list);
            catch exception
                bad = find(or(state_n_list > length(obj.ValueMap), state_n_list < 1), 1);
                error('rpMarkov:InvalidStateNum', ...
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
                error('rpMarkov:InvalidValue', ...
                    'Attempt to find value for invalid value: %g', ...
                    values(bad))
            end
            if size(values) ~= size(state_n_list)
                state_n_list = reshape(state_n_list, size(values));
            end
        end

        function [value_series, state_n_series] = dsim(obj, t_list, initial_value)
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
        % Note: after calling dsim, the process internal time will be set to
        % the final value of t_list

            %identify valid simulation times
            ok = (t_list >= 0);

            %initialize outputs
            value_series = zeros(size(t_list));
            value_series(not(ok)) = NaN;
            state_n_series = zeros(size(value_series));
            state_n_series(not(ok)) = NaN;

            %only simulate valid values of t_list
            t_list = t_list(ok);

            %only run the simulation if there are valid times to simulate.
            %If not, we have already filled the value list with NaNs and
            %can skip ahead to the state list if requested.
            if not(isempty(t_list))
                %Find time range that we need to simulate
                t_max = max(t_list);
                t_min = min(t_list);

                %if we are given an initial value, check that it is valid
                %for the given minimum simulation time and start there
                if nargin >= 3 && not(isempty(initial_value))
                    cur_idx = obj.s2idx(obj.dval2num(initial_value),t_min);
                    if cur_idx == 0
                        error('rpMarkov:InvalidValueAtTime', ...
                            'Initial value, %g, not valid at first time, %d', ...
                            initial_value, t_min)
                    end
                    obj.cur_state = obj.LatticeSNum{t_min+1}(cur_idx);
                else
                    %if all simulation times are after our current time, use the
                    %current state and time as a starting point, but if not, reset
                    %to the starting point and time.
                    if obj.t > t_min
                        obj.reset()
                        t_min = 0;
                    else
                        t_min = floor(obj.t);
                    end
                end

                %initialize output vectors
                t_sim = t_min:t_max;
                %make sure we match the required output vector shape
                if size(t_list,2) == 1
                    t_sim = t_sim';
                end

                v_sim = zeros(size(t_sim));
                s_sim = zeros(size(t_sim));

                %store initial time step
                v_sim(1) = obj.dnum2val(obj.cur_state);
                s_sim(1) = obj.cur_state;

                %simulate future timesteps
                for idx = 2:length(t_sim)
                    [s_sim(idx), v_sim(idx)] = obj.step_helper(s_sim(idx-1), t_sim(idx));
                end

                %Stuff the correct values into the output list
                value_series(ok) = v_sim(t_list-t_min+1);
                state_n_series(ok) = s_sim(t_list-t_min+1);

                obj.t = t_list(end);
                obj.cur_state = s_sim(obj.t-t_min+1);
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

            if ischar(t) && strcmp(t, 'all')
                value_range = obj.ValueMap([1, end]);
                state_n_range = [1, length(obj.ValueMap)];
            else
                v_idx = obj.t2idx(t);

                value_range = [min(obj.VSet{v_idx}), max(obj.VSet{v_idx})];
                if nargout > 1
                    % Note: this formulation relies on the fact that the
                    % SSet and VSet are both in the same order
                    state_n_range = [min(obj.SSet{v_idx}), max(obj.SSet{v_idx})];
                end
            end
        end


        %% ===== Additional simulation support
        function [value, state_n, t] = step(obj, delta_t)
        %STEP simulate forward or backward
        %
        % by default steps forward by delta_t = 1
            if nargin < 2
                delta_t = 1;
            end
            %compute the proposed new time
            new_t = obj.t + delta_t;

            %check if it is valid, if not return empty results
            if new_t < 0
                value = [];
                state_n = [];
                t = obj.t;
                return
            else
                %if new time is valid simulate forward or back as needed
                if floor(obj.t) == floor(new_t)
                    %No need to actually change, b/c we have discrete steps
                    %that round
                    value = obj.dnum2val(obj.cur_state);
                else
                    if new_t > obj.t
                        %Step forward
                        for t = floor(obj.t):(floor(new_t)-1)
                            obj.cur_state = obj.step_helper(obj.cur_state, t);
                        end
                    else
                        %Step backward
                        for t = floor(obj.t):-1:(floor(new_t)+1)
                            %Extract possible previous states (discrete)
                            [~, state_n_list, prob] = obj.dlistprev(obj.cur_state, t );

                            % Randomly select the new state based on the cdf
                            % described by the prob vector (computed using cumsum)
                            new_idx = find(rand <= cumsum(prob), 1, 'first');
                            obj.cur_state = state_n_list(new_idx);
                        end
                    end

                    % Update our stored state
                    obj.t = t;

                    % And finally setup our output
                    state_n = obj.cur_state;

                end

            end

        end

        function [value, state_n, t] = curState(obj)
        %CURSTATE Return the current state of the simulation
            value = obj.ValueMap(obj.cur_state);
            state_n = obj.cur_state;
            t = obj.t;
        end

        function reset(obj, initial_value)
        %RESET reset simulation to t=0 and provided initial_value. If no
        %initial value provided, we randomly sample the first state
            obj.t = 0;
            if nargin > 1
                obj.cur_state = obj.dval2num(initial_value);
            else
                obj.cur_state = obj.SSet{1}(find(rand() <= obj.CdfSet{1}, 1, 'first'));
            end
        end

    end

end
