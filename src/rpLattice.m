classdef rpLattice < RandProcess
%rpLattice bi/multi-nomial lattice random process: Geometric Brownian Motion
%
% Random process that models a discrete random walk with constant geometric
% factors such that the branches of the tree recombine such that each time
% step's outcome space grows linearly rather than exponentially.
%
% Constructor: obj = rpLattice(start, coef, prob, t_max, tol)
%
% Required inputs/properties:
%  Start:  starting value
%  Coef: vector of growth coeficients
%  Prob:   vector of corresponding coefficients
%  MaxT: maximum time for lattice (starting from t=0)
% Optional inputs/properties:
%  Tol:  rounding tolerance for combining states (default = 0.0001)
%
% Notes:
%  - Integer timesteps are assumed, scale accordingly for fractional
%    timesteps
%  - start is the value at t=0
%  - For non-integer times, the output is assumed to remain at the value
%    corresponding to the previous integer time step (Zero-order hold)
%  - For t>MaxT, the lattice is assumed to remain constant with no more
%    transitions
%
% see also RandProc, rpDiscreteSample, rpMarkov, rpList
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-12-16 14:30  BryanP      Adapted from rpList v6
%   2  2010-12-15 23:30  BryanP      Added t_max for dlist(... 'all')
%   3  2010-12-21 21:20  BryanP      Made Tmax required in constructor, etc
%   4  2010-12-23 15:45  BryanP      Implemented buildLattice & dlist*
%   5  2011-01-04 22:00  BryanP      Complete, working version
%   6  2011-03-12 12:45  BryanP      Reordered set.Coef & Prob to work with save/load
%   7  2011-06-09 16:00  BryanP      match input shape for dnum2val and dval2num
%   8  2011-11-09 23:40  BryanP      Added tolerance to constructor
%   9  2012-01-11 17:00  BryanP      Added constant non-walk for t>Tmax
%  10  2012-01-25 13:45  BryanP      Allow blank constructor input for loading from a file
%  11  2012-04-17 08:55  BryanP      Reworked cdf for sim/sample

    properties
        Start = [];   % initial value
        Coef = [];    % vector of coefficients
        Prob = [];    % vector of probabilities
        Tol = 0.0001; % rounding tolerance for merging states
        Tmax = 0;     % maximum time for which to compute the lattice
        t = 0;        % Current time
    end

    % Internal properties
    properties (Access='protected')
        LatticeValue = {}  %Cell array of precomputed values for each t
        LatticeProb  = {}  %Cell array of precomuted probabilities at t
        LatticeSNum = {}   %Cell array of precomputed state numbers for each t

        ValueMap = []    %List of unique values in lattice, index provides state number

        cur_state = NaN    %Current state number for simulation

        %Flags
        SkipLatticeBuild = false    %Wait till all values set before building the lattice
        SuppressLengthError = false %Prevent coef & prob length checks when changing both

    end

    methods (Access = protected)
        %Reset the stored Lattice
        function clearStoredLattice(obj)
            obj.LatticeValue = {obj.Start};
            obj.LatticeSNum = {1};
            obj.LatticeProb = {1};
            obj.ValueMap = 1;
        end

        %Build & cache the lattice
        function buildLattice(obj)
            if obj.SkipLatticeBuild
                return
            end
            % -- Build the value lattice.
            %Starting with the initial value
            obj.LatticeValue = {obj.Start};
            obj.LatticeProb = {1};
            %Then cycling through future time scenarios
            for idx=2:(obj.Tmax+1)
                % Compute next lattice state by
                % 1) multiplying the previous set of states by the
                % coefficient matrix to form an array
                %Implementation Note: have to put coef first, otherwise,
                %the first multiplication creates a row vector and unique
                %does not change it back to a column vector.
                obj.LatticeValue{idx} = obj.Coef * obj.LatticeValue{idx-1}';
                obj.LatticeProb{idx} = obj.Prob * obj.LatticeProb{idx-1}';

                % 2) rounding these values to make sure that states merge
                % Note: the RoundTo function is part of the adp toolbox
                obj.LatticeValue{idx} = RoundTo(obj.LatticeValue{idx}, obj.Tol);

                % 3) converting the resulting rectangular matrices to
                % column vectors.
                %
                % Note: in reshape, the [] will compute the required number
                % of rows.
                obj.LatticeValue{idx} = reshape(obj.LatticeValue{idx},[],1);
                obj.LatticeProb{idx} = reshape(obj.LatticeProb{idx},[],1);

                % 4) Using unique() to return an ordered vector of unique
                % values.
                [obj.LatticeValue{idx}, junk, idx_list]  = unique(obj.LatticeValue{idx}); %#ok<ASGLU>

                % 5) Use the indices (idx_list) to locate the corresponding
                % probabilities and then sum the combined probabilities
                % using accumarray
                obj.LatticeProb{idx} = accumarray(idx_list, obj.LatticeProb{idx});

                % 6) Normalize the probability so we are sure to sum to 1
                obj.LatticeProb{idx} = obj.LatticeProb{idx}/sum(obj.LatticeProb{idx});

            end

            % -- Now compute the list of possible values
            % Note: the corresponding state number is simply its vector
            % index
            obj.ValueMap = unique(vertcat(obj.LatticeValue{:}));

            % -- And store the state list matrix
            for idx=1:(obj.Tmax+1)
                obj.LatticeSNum{idx} = obj.dval2num(obj.LatticeValue{idx});
            end

        end
    end

    methods
        %% ===== Constructor & related
        function obj = rpLattice(start, coef, prob, t_max, tol)
            % Support zero parameter calls to constructor for special
            % MATLAB situations (see help files)
            if nargin > 0
                switch nargin
                    case 5
                        obj.setparams(start, coef, prob, t_max, tol);
                    case 4
                        obj.setparams(start, coef, prob, t_max);
                    otherwise
                        obj.setparams(start, coef, prob);
                end
            end
        end

        function setparams(obj, start, coef, prob, t_max, tol)
        %SETPARAMS allow easy changing of all lattice parameters at once
            %Prevent building the stored lattice until required values set
            obj.SkipLatticeBuild = true;
            %Prevent errors from set functions for length mis-matches
            obj.SuppressLengthError = true;

            if nargin > 4
                obj.Tmax = t_max;
            end
            if nargin > 5
                obj.Tol = tol;
            end
            obj.Start = start;
            %We have checks to make sure that coef & prob are the same
            %length, since we have both supressed
            if length(coef) ~= length(prob)
                error('rpLattice:CoefProbMismatch', 'Coefficient & Probability vectors must have equal length')
            end

            %Now set probability
            obj.Prob = prob;

            %Let set.Coef handle precomputing the lattice & return lattice
            %build state to default
            obj.SkipLatticeBuild = false;

            %And finally set coefficient list
            obj.Coef = coef;

            % Reenable coef & prob length checking
            obj.SuppressLengthError = false;

            obj.reset()
        end


        %% ===== Property value maintenance
        % Clear stored values & probs when changing start value
        % Also check that the start value is scalar
        function set.Start(obj, s)
            if isscalar(s)
                obj.Start = s;
                obj.buildLattice();
            else
                %Allow empty for loading Lattice objects from *.mat files
                if not(isempty(s))
                    error('rpLattice:NonScalarStart', 'Lattice starting value must be a scalar')
                end
            end
        end

        % Clear stored values & probs when changing coeficient vector
        % Also check that coef & prob are the same length
        function set.Coef(obj, c)
            % Check for only positive (or zero) coefficient values
            if any(c<0)
                error('rpLattice:NegCoef', 'All coeficients must be positive')
            end

            %convert to column vector if needed
            if size(c, 2) ~= 1;
                c = c';
            end
            obj.Coef = c;

            % Ensure coef length matches prob length, unless this check is
            % suppressed b/c setting both
            %
            % Note: we do this after assigning obj.Prob to enable
            % successful saves & loads.
            if not(obj.SuppressLengthError) %#ok<MCSUP>
                if length(c) ~= length(obj.Prob) %#ok<MCSUP>
                        %Allow empty for loading Lattice objects from *.mat files
                        if not(isempty(c)) && not(isempty(obj.Prob))
                            error('rpLattice:CoefProbMismatch', ...
                                'Coefficient and Probability vectors must have equal length')
                        end
                end
            end

            %Rebuild the stored lattice
            obj.buildLattice;
        end

        % Clear stored values & probs when changing probability vector
        % Also check that coef & prob are the same length
        function set.Prob(obj, p)
            %Check that we have a valid probability vector
            if sum(p) ~= 1
                %Allow empty for loading Lattice objects from *.mat files
                if not(isempty(p))
                    error('rpLattice:ProbNotSumOne', 'Probability must sum to one')
                end
            end

            %convert to column vector if needed
            if size(p, 2) ~= 1;
                p = p';
            end
            obj.Prob = p;

            % Ensure prob length matches coef length, unless this check is
            % suppressed b/c setting both
            %
            % Note: we do this after assigning obj.Prob to enable
            % successful saves & loads.
            if not(obj.SuppressLengthError) %#ok<MCSUP>
                if length(p) ~= length(obj.Coef) %#ok<MCSUP>
                    error('rpLattice:CoefProbMismatch', ...
                        'Probability and Coefficient vectors must have equal length')
                end
            end

            %Rebuild lattice
            obj.buildLattice;
        end

        % If changing t_max, remove no longer needed stored Lattice data
        function set.Tmax(obj, t_max)
            obj.Tmax = t_max;
            obj.buildLattice;
        end

        % Check that we set a valid t
        function set.t(obj, t)
            if t < 0;
                error('RandProcess:InvalidTime', 'time (t) must be positive')
            else
                if obj.t ~= t
                    obj.t = t;

                    %Now set the current state to the middle state for this
                    %time and limit the value to Tmax
                    t = min(floor(t),obj.Tmax); %#ok<MCSUP>
                    s_idx = ceil(length(obj.LatticeSNum{t+1})/2); %#ok<MCSUP>
                    obj.cur_state = obj.LatticeSNum{t+1}(s_idx); %#ok<MCSUP>
                end
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
        % case, t_max is REQUIRED to bound the problem size.
            if nargin < 2 || isempty(t)
                [value_list, state_n_list] = obj.dlist(obj.t);
            elseif (ischar(t) && strcmp(t, 'all'))
                value_list = obj.ValueMap;
                state_n_list = 1:length(obj.ValueMap);
            elseif t < 0;
                error('RandProcess:InvalidTime', 'Only t>0 valid for rpLattice')
            else
                t = min(obj.Tmax, t);
                value_list = obj.LatticeValue{t+1};
                state_n_list = obj.dval2num(value_list);
            end
        end

        function [value_list, state_n_list, prob] = dlistprev (obj, state_n, t )
        % DLISTPREV List previous discrete states & probabilities
        %
        % List possible previous states (by number) along with conditional
        % probability P(s_{t-1} | s_t)
        %
        % If t is not defined, the current simulation time is assumed
            if nargin < 3
                if nargin < 2
                    state_n = obj.cur_state;
                end
                [value_list, state_n_list, prob] = obj.dlistprev(state_n, obj.t);
                return
            elseif t <= 0;
                error('RandProcess:InvalidTime', 'Only t>0 valid for rpLattice')
            end

            %find a valid time for state lookup, by limiting to Tmax
            t_lookup = min(t, obj.Tmax);

            if isempty(state_n) ...
                    || state_n > length(obj.ValueMap) ...
                    || not(all(ismember(obj.ValueMap(state_n),obj.LatticeValue{t_lookup+1})))
                error('RandProcess:InvalidState', 'State #%d not valid at time=%d', state_n, t)
            else
                %if we get here, we know the state is valid for this time

                if t > obj.Tmax
                    value_list = obj.ValueMap(state_n);
                else

                    % Build previous value list by reversing the multiplication
                    % by coef required to get there
                    value_list = RoundTo(obj.ValueMap(state_n) ./ obj.Coef, obj.Tol);

                    % For values near the "edge" not all of the possible priors
                    % from this division are actually valid states during the
                    % previous period, so limit our seach to those that are.
                    % Note: that the index t, corresponds to t-1 b/c 1-indexed
                    [value_list, coef_used, states] = intersect(value_list, obj.LatticeValue{t});
                end

                if nargout >1
                    %Now we have enough information to return the state list
                    state_n_list = obj.dval2num(value_list);
                end

                if nargout > 2
                    if t > obj.Tmax
                        prob = 1;
                    else
                        %Compute the probability using Bayes' Theorem:
                        %                     P(s_t | s_{t-1}) P(s_{t-1})
                        % P(s_{t-1} | s_t) =  ---------------------------
                        %                           P{s_t)
                        %
                        % But rather than explicitly computing P{s_t}, simply
                        % normalize by dividing by the sum of the probabilities
                        % which is equivalent.
                        prob = obj.Prob(coef_used) .* obj.LatticeProb{t}(states);
                        prob = prob/sum(prob);
                    end
                end
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
                [value_list, state_n_list, prob] = obj.dlistnext(state_n, obj.t);
                return
            elseif t < 0;
                error('RandProcess:InvalidTime', 'Only t>0 valid for rpLattice')
            end

            %find a valid time for state lookup, by limiting to Tmax
            t_lookup = min(t, obj.Tmax);

            if isempty(state_n) ...
                    || state_n > length(obj.ValueMap) ...
                    || not(all(ismember(obj.ValueMap(state_n),obj.LatticeValue{t_lookup+1})))
                error('RandProcess:InvalidState', 'State #%d not valid at time=%d', state_n, t)
            else
                %if we get here, we know the state is valid for this time

                if t > obj.Tmax -1
                    value_list = obj.ValueMap(state_n);
                else
                    % Build next value list by multiplying
                    % by coef required to get there
                    value_list = RoundTo(obj.ValueMap(state_n) .* obj.Coef, obj.Tol);
                end

                if nargout >1
                    %Now we have enough information to return the state list
                    state_n_list = obj.dval2num(value_list);
                end

                if nargout > 2
                    %Here the probability is easy... it is either
                    if t > obj.Tmax -1
                        %one or
                        prob = 1;
                    else
                        %our probability vector
                        prob = obj.Prob;
                    end
                end
            end
        end

        function values = dnum2val (obj, state_n_list )
        % DNUM2VAL Convert discrete state number to a value
            try
                values = obj.ValueMap(state_n_list);
            catch exception
                bad = find(or(state_n_list > length(obj.ValueMap), state_n_list < 1), 1);
                error('rpLattice:InvalidStateNum', ...
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
                error('rpLattice:InvalidValue', ...
                    'Attempt to find state number for invalid value: %g', ...
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
        % Invalid times (t<0 or t>Tmax) return NaN
        %
        % Note: after calling dsim, the process internal time will be set to
        % the final value of t_list

            %identify simulation time classes
            ok = (t_list >= 0);
            above_Tmax = t_list > obj.Tmax;
            ok_below_Tmax = ok & not(above_Tmax);

            %initialize outputs
            value_series = zeros(size(t_list));
            value_series(not(ok)) = NaN;

            %only simulate valid values of t_list
            t_ok = t_list(ok);
            t_ok_below_Tmax = t_list(ok_below_Tmax);

            %only run the simulation if there are valid times to simulate.
            %If not, we have already filled the value list with NaNs and
            %can skip ahead to the state list if requested.
            if not(isempty(t_ok))
                %Find time range that we need to simulate
                t_max = min(obj.Tmax, max(t_ok));
                t_min = min(t_ok);

                %if we are given an initial value, check that it is valid
                %for the given minimum simulation time and start there
                if nargin >= 3 && not(isempty(initial_value))
                    cur_idx = find(initial_value == obj.LatticeValue{t_min+1}, 1);
                    if isempty(cur_idx)
                        error('rpLattice:InvalidValueAtTime', ...
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

                %compute transition distribution function.
                trans_cdf = cumsum(obj.Prob);

                %store initial time step
                v_sim(1) = obj.dnum2val(obj.cur_state);
                %simulate future timesteps
                for idx = 2:length(t_sim)
                    trans = find(rand() <= trans_cdf, 1, 'first');
                    v_sim(idx) = obj.Coef(trans) * v_sim(idx-1);
                    v_sim(idx) = RoundTo(v_sim(idx),obj.Tol);
                end

                %Stuff the correct values into the output list
                value_series(ok_below_Tmax) = v_sim(t_ok_below_Tmax-t_min+1);
                value_series(above_Tmax) = v_sim(end);

                obj.t = t_ok(end);
                value_series_ok = value_series(ok);
                obj.cur_state = obj.dval2num(value_series_ok(end));
            end

            %Handle state numbers if required
            if nargout > 1;
                state_n_series = zeros(size(value_series));
                state_n_series(not(ok)) = NaN;
                state_n_series(ok) = obj.dval2num(value_series(ok));
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
                value_range = [obj.ValueMap(1), obj.ValueMap(end)];
                state_n_range = [1, length(obj.ValueMap)];
            else
                %Handle any non integer values
                t = floor(t);

                if t < 0;
                    error('RandProcess:InvalidTime', 'Only t>0 valid for rpLattice')
                end

                t = min(t, obj.Tmax);
                value_range = [min(obj.LatticeValue{t+1}), max(obj.LatticeValue{t+1})];
                if nargout > 1
                    state_n_range = obj.dval2num(value_range);
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

            %check if it is valid
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
                    state_n = obj.cur_state;
                    value = obj.dnum2val(state_n);
                else
                    if new_t > obj.t
                        %Step forward
                        trans_prob = cumsum(obj.Prob);
                        for t = floor(obj.t):(floor(new_t)-1)
                            if t >= obj.Tmax
                                value_list = obj.ValueMap(obj.cur_state);
                                new_idx = 1;
                            else
                                %Extract possible next states (discrete)
                                value_list = RoundTo(obj.ValueMap(obj.cur_state) .* obj.Coef, obj.Tol);

                                % Randomly select the new state based on the cdf
                                % described by the prob vector (computed using cumsum)
                                new_idx = find(rand<=trans_prob, 1, 'first');
                                obj.cur_state = obj.dval2num(value_list(new_idx));
                            end
                        end
                    else
                        %Step backward
                        for t = floor(obj.t):-1:(floor(new_t)+1)
                            %Extract possible previous states (discrete)
                            [value_list, state_n_list, prob] = obj.dlistprev(obj.cur_state, t );

                            % Randomly select the new state based on the cdf
                            % described by the prob vector (computed using cumsum)
                            new_idx = find(rand <= cumsum(prob), 1, 'first');
                            obj.cur_state = state_n_list(new_idx);
                        end
                    end
                    value = value_list(new_idx);
                    state_n = obj.cur_state;
                end

            end
            % Ensure that we will always leave the object's time in a valid
            % state by truncating excessively high or small stepsizes
            obj.t = max(0, new_t);
            t = obj.t;

            % Setting t changes to the default state for the current time
            % so make sure we leave ourselves in the actual simulated
            % state.
            obj.cur_state = state_n;
        end

        function [value, state_n, t] = curState(obj)
        %CURSTATE Return the current state of the simulation
            value = obj.ValueMap(obj.cur_state);
            state_n = obj.cur_state;
            t = obj.t;
        end

        function reset(obj, initial_value) %#ok<INUSD>
        %RESET reset simulation to t=0
            obj.cur_state = obj.LatticeSNum{1};
            obj.t = 0;

        end

    end

end
