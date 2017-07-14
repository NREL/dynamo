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
%  MaxT: maximum time for lattice (starting from t=1)
% Optional inputs/properties:
%  Tol:  rounding tolerance for combining states (default = 0.0001)
%
% Notes:
%  - Integer timesteps are assumed, scale accordingly for fractional
%    timesteps
%  - start is the value at t=1
%  - For non-integer times, the output is assumed to remain at the value
%    corresponding to the previous integer time step (Zero-order hold)
%  - For t>MaxT, the lattice is assumed to remain constant with no more
%    transitions
%
% see also rpPathDepend, RandProc, rpDiscreteSample, rpMarkov, rpBasic
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  13  2017-07-14 05:48  BryanP      Only support value states 
%  12  2017-07-13 21:17  BryanP      Initial Update for new RandProc format
%  11  2012-04-17 08:55  BryanP      Reworked cdf for sim/sample
%  10  2012-01-25 13:45  BryanP      Allow blank constructor input for loading from a file
%   9  2012-01-11 17:00  BryanP      Added constant non-walk for t>Tmax
%   8  2011-11-09 23:40  BryanP      Added tolerance to constructor
%   7  2011-06-09 16:00  BryanP      match input shape for dnum2val and dval2num
%   6  2011-03-12 12:45  BryanP      Reordered set.Coef & Prob to work with save/load
%   5  2011-01-04 22:00  BryanP      Complete, working version
%   4  2010-12-23 15:45  BryanP      Implemented buildLattice & dlist*
%   3  2010-12-21 21:20  BryanP      Made Tmax required in constructor, etc
%   2  2010-12-15 23:30  BryanP      Added t_max for dlist(... 'all')
%   1  2010-12-16 14:30  BryanP      Adapted from rpList v6

    properties
        Start = [];   % initial value
        Coef = [];    % vector of coefficients
        Prob = [];    % vector of probabilities
        Tol = 0.0001; % rounding tolerance for merging states
        Tmax = 0;     % maximum time for which to compute the lattice
        t = 1;        %current timestep

    end

    % Internal properties
    properties (Access='protected')
        LatticeValue = {}  %Cell array of precomputed values for each t
        LatticeUncondProb  = {}  %Cell array of precomuted (unconditional) probabilities at t
        LatticeUncondCdf  = {}   %Cell array of precomuted (unconditional) cdfs at t
        ConditionalCdf = [];          %precomputed conditional probability (based on Prob)
        
        %Flags
        SkipLatticeBuild = false    %Wait till all values set before building the lattice
        SuppressLengthError = false %Prevent coef & prob length checks when changing both
    end

    methods (Access = protected)
        %Reset the stored Lattice
        function clearStoredLattice(obj)
            obj.LatticeValue = {obj.Start};
            obj.LatticeUncondProb = {1};
            obj.LatticeUncondCdf = {1};
        end

        %Build & cache the lattice
        function buildLattice(obj)
            if obj.SkipLatticeBuild
                return
            end
            % -- Build the value lattice.
            %compute conditional cdf here to avoid parameter
            %cross-reference issues
            obj.ConditionalCdf = cumsum(obj.Prob);

            %Starting with the initial value
            obj.LatticeValue = {obj.Start};
            obj.LatticeUncondProb = {1};
            obj.LatticeUncondCdf = {1};
            %Then cycling through future time scenarios
            for idx=2:(obj.Tmax+1)
                % Compute next lattice state by
                % 1) multiplying the previous set of states by the
                % coefficient matrix to form an array & compute
                % corresponding unconditional probability (for sampling)
                %Implementation Note: have to put coef first, otherwise,
                %the first multiplication creates a row vector and unique
                %does not change it back to a column vector.
                obj.LatticeValue{idx} = obj.Coef * obj.LatticeValue{idx-1}';
                obj.LatticeUncondProb{idx} = obj.Prob * obj.LatticeUncondProb{idx-1}';

                % 2) rounding these values to make sure that states merge
                % Note: the RoundTo function is part of the adp toolbox
                obj.LatticeValue{idx} = RoundTo(obj.LatticeValue{idx}, obj.Tol);

                % 3) converting the resulting rectangular matrices to
                % column vectors.
                %
                % Note: in reshape, the [] will compute the required number
                % of rows.
                obj.LatticeValue{idx} = reshape(obj.LatticeValue{idx},[],1);
                obj.LatticeUncondProb{idx} = reshape(obj.LatticeUncondProb{idx},[],1);

                % 4) Using unique() to return an ordered vector of unique
                % values.
                [obj.LatticeValue{idx}, junk, idx_list]  = unique(obj.LatticeValue{idx}); %#ok<ASGLU>

                % 5) Use the indices (idx_list) to locate the corresponding
                % probabilities and then sum the combined probabilities
                % using accumarray
                obj.LatticeUncondProb{idx} = accumarray(idx_list, obj.LatticeUncondProb{idx});

                % 6) Normalize the probability so we are sure to sum to 1
                obj.LatticeUncondProb{idx} = obj.LatticeUncondProb{idx}/sum(obj.LatticeUncondProb{idx});
                
                % 7) Compute corresponding unconditional cdf
                obj.LatticeUncondCdf{idx} = cumsum(obj.LatticeUncondProb{idx});
            end
        end
    end

    methods
        %% ===== Constructor & related
        function obj = rpLattice(start, coef, prob, t_max, tol)
            % Note: see rpLattice class documentation (above) for more info
            %
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
            obj.reset();
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
            if isrow(c)
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
                        if not(isempty(c)) && not(isempty(obj.Prob)) %#ok<MCSUP> OK b/c we explicitly check for empty (unset) values
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
            if isrow(p)
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
                    %Allow empty for loading Lattice objects from *.mat files
                    if not(isempty(p)) && not(isempty(obj.Coef)) %#ok<MCSUP> OK b/c we explicitly check for empty (unset) values
                        error('rpLattice:CoefProbMismatch', ...
                            'Coefficient and Probability vectors must have equal length')
                    end
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
            if t < 1
                error('RandProcess:InvalidTime', 'time (t) must be strictly positive')
            else
                if obj.t ~= t
                    obj.t = t;

                    %Now set the current state to the middle state for this
                    %time and limit the value to Tmax
                    t = min(floor(t),obj.Tmax); %#ok<MCSUP>
                    s_idx = ceil(length(obj.LatticeSNum{t+1})/2);
                    obj.cur_state = obj.LatticeSNum{t+1}(s_idx);
                end
            end
        end


        %% == Additional public methods
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
        %   state_list = sample(obj, N, t, cur_state)
        %       Sample specified time using conditional probability
        %       starting from provided state
 
            if nargin < 2
                N = 1;
            end
            
            if nargin < 3 || isempty(t)
                t=obj.t;
            end
            
            if nargin < 4 || isempty(cur_state)
                cur_state=obj.cur_state;
            end
            
            %Handle any non integer or large values for time
            t = min(floor(t), obj.Tmax);

            idx_at_t = zeros(N,1);
            for samp_idx = 1:N
                idx_at_t(samp_idx) = find(rand(1) <= obj.CdfList{t}, 1, 'first');
            end
            state_list = obj.Vlist{t}(idx_at_t,:);
        end
        
        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.

        function state_list = dlist (obj, t)
        % DLIST List possible discrete states
        %
        % List possible discrete states by number for given time
        % if t is not listed, the states for the current simulation time
        % are returned.
        %
        % To get a list of all possible states pass with t='all'
            if nargin < 2 || isempty(t)
                state_list = obj.dlist(obj.t);
            elseif (ischar(t) && strcmp(t, 'all'))
                state_list = obj.ValueMap;
            elseif t < 1
                error('RandProcess:InvalidTime', 'Only t>1 valid for rpLattice')
            else
                t = min(obj.Tmax, t);
                state_list = obj.LatticeValue{t+1};
            end
        end

        function [state_list, prob] = dlistprev (obj, state_in, t )
        % DLISTPREV List previous discrete states & probabilities
        %
        % List possible previous states along with conditional
        % probability P(s_{t-1} | s_t)
        %
        % If t is not defined, the current simulation time is assumed
            if nargin < 3
                if nargin < 2
                    state_in = obj.cur_state;
                end
                [state_list, prob] = obj.dlistprev(state_in, obj.t);
                return
            elseif t <= 1
                error('RandProcess:InvalidTime', 'Only t>1 valid for rpLattice')
            end

            %find a valid time for state lookup, by limiting to Tmax
            t_lookup = min(t, obj.Tmax);

            if isempty(state_in) ...
                    || not(all(ismembertol(state_in,obj.LatticeValue{t_lookup+1}, obj.Tol)))
                error('RandProcess:InvalidState', 'State %g not valid at time=%d', state_in, t)
            else
                %if we get here, we know the state is valid for this time

                if t > obj.Tmax
                    state_list = state_in;
                else

                    % Build previous value list by reversing the multiplication
                    % by coef required to get there
                    state_list = RoundTo(obj.ValueMap(state_in) ./ obj.Coef, obj.Tol);

                    % For values near the "edge" not all of the possible priors
                    % from this division are actually valid states during the
                    % previous period, so limit our seach to those that are.
                    % Note: that the index t, corresponds to t-1 b/c 1-indexed
                    [state_list, coef_used, states] = intersect(state_list, obj.LatticeValue{t});
                end

                if nargout > 1
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
                        prob = obj.Prob(coef_used) .* obj.LatticeUncondProb{t}(states);
                        prob = prob/sum(prob);
                    end
                end
            end
        end

        function [state_list, prob] = dlistnext (obj, state_in, t )
        % DLISTNEXT List next discrete states & probabilities
        %
        % List possible next states (by number) along with conditional
        % probability P(s_{t+1} | s_t)
        %
        % If t and/or state are not defined, the current simulation time
        % and state are assumed
            if nargin < 3
                [state_list, prob] = obj.dlistnext(state_in, obj.t);
                return
            elseif t < 1
                error('RandProcess:InvalidTime', 'Only t>0 valid for rpLattice')
            end

            %find a valid time for state lookup, by limiting to Tmax
            t_lookup = min(t, obj.Tmax);

            if isempty(state_in) ...
                    || not(all(ismembertol(state_in,obj.LatticeValue{t_lookup+1}, obj.Tol)))
                error('RandProcess:InvalidState', 'State %g not valid at time=%d', state_in, t)
            else
                %if we get here, we know the state is valid for this time

                if t > obj.Tmax
                    state_list = obj.ValueMap(state_in);
                else
                    % Build next value list by multiplying
                    % by coef required to get there
                    state_list = RoundTo(state_in .* obj.Coef, obj.Tol);
                end

                if nargout > 2
                    %Here the probability is easy... it is either
                    if t > obj.Tmax
                        %one or
                        prob = 1;
                    else
                        %our probability vector
                        prob = obj.Prob;
                    end
                end
            end
        end

        function value_series = dsim(obj, t_list, initial_value)
        % DSIM Simulate discrete process.
        %
        % A column vector for t is assumed to be a series of times for
        % which to return results. Intermediate times are also computed, if
        % needed, but not returned. The initial value is not returned in
        % the value series. Only one simulation is run, such that out of
        % order times will be sorted before simulation and duplicate times
        % will return the same result
        %
        % Invalid times (t<1 or t>Tmax) return NaN
        %
        % Note: after calling dsim, the process internal time will be set to
        % the final value of t_list

            %identify simulation time classes
            ok = (t_list > 0);
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
                    initial_value = RoundTo(initial_value, obj.Tol);
                    cur_idx = find(initial_value == obj.LatticeValue{t_min}, 1);
                    if isempty(cur_idx)
                        error('rpLattice:InvalidValueAtTime', ...
                            'Initial value, %g, not valid at first time, %d', ...
                            initial_value, t_min)
                    end
                    obj.cur_state = initial_value;
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
                if isrow(t_list)
                    t_sim = t_sim';
                end

                v_sim = zeros(size(t_sim));

                %compute transition distribution function.
                trans_cdf = cumsum(obj.Prob);

                %store initial time step
                v_sim(1) = obj.cur_state;
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
                obj.cur_state = value_series_ok(end);
            end

        end

        %% ===== General (discrete or continuous) Methods
        function value_series = sim(obj, t_list, initial_value)
        % SIM Simulate process for desired (continuous) times
        %
        % A column vector for t is assumed to be a series of times for
        % which to return results. Intermediate times are also computed, if
        % needed. The initial value is not returned in the value series.
        %
        % Function must handle arbitrary positive values for t_list
        % Invalid times (t<=1) return NaN
        %
        % Note: after calling sim, the process internal time will be set to
        % the final value of t_list
        %
        % This implementation typically rounds down to the nearest integer
        % (zero order hold)
            if nargin < 3
                initial_value = [];
            end
            value_series = obj.dsim(floor(t_list), initial_value);

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
                value_range = [min(vertcat(obj.LatticeValue{:})), max(vertcat(obj.LatticeValue{:}))];
            else
                %Handle any non integer values
                t = floor(t);

                if t < 1
                    error('RandProcess:InvalidTime', 'Only t>=1 valid for rpLattice')
                end

                t = min(t, obj.Tmax);
                value_range = [min(obj.LatticeValue{t}), max(obj.LatticeValue{t})];
            end
        end


        %% ===== Additional simulation support
        function [value, t] = step(obj, delta_t)
        %STEP simulate forward or backward
        %
        % by default steps forward by delta_t = 1
            if nargin < 2
                delta_t = 1;
            end
            %compute the proposed new time
            new_t = obj.t + delta_t;

            %check if it is valid
            if new_t < 1
                value = [];
                t = obj.t;
                return
            else
                %if new time is valid simulate forward or back as needed
                if floor(obj.t) == floor(new_t)
                    %No need to actually change, b/c we have discrete steps
                    %that round
                    value = obj.cur_state;
                else
                    if new_t > obj.t
                        %Step forward
                        for t = floor(obj.t):(floor(new_t))
                            if t >= obj.Tmax
                                state_list = obj.cur_state;
                                new_idx = 1;
                            else
                                trans_prob = cumsum(obj.LatticeUncondProb{t});
                                %Extract possible next states (discrete)
                                state_list = RoundTo(obj.cur_state .* obj.Coef, obj.Tol);

                                % Randomly select the new state based on the cdf
                                % described by the prob vector (computed using cumsum)
                                new_idx = find(rand<=trans_prob, 1, 'first');
                                obj.cur_state = state_list(new_idx);
                            end
                        end
                    else
                        %Step backward
                        for t = floor(obj.t):-1:(floor(new_t))
                            %Extract possible previous states (discrete)
                            [state_list, prob] = obj.dlistprev(obj.cur_state, t );

                            % Randomly select the new state based on the cdf
                            % described by the prob vector (computed using cumsum)
                            new_idx = find(rand <= cumsum(prob), 1, 'first');
                            obj.cur_state = state_list(new_idx);
                        end
                    end
                    value = state_list(new_idx);
                end

            end
            % Ensure that we will always leave the object's time in a valid
            % state by truncating excessively high or small stepsizes
            obj.t = max(1, new_t);
            t = obj.t;

            % Setting t changes to the default state for the current time
            % so make sure we leave ourselves in the actual simulated
            % state.
            obj.cur_state = value;
        end

        function [value, t] = cur_state(obj)
        %CURSTATE Return the current state of the simulation
            value = obj.cur_state;
            t = obj.t;
        end

        function reset(obj, initial_value) %#ok<INUSD>
        %RESET reset simulation to t=1
            if nargin > 1
                warning('rpLattice:InvalidStateForTime', 'Cannot specify value for rpLattice reset (only one t=1 value possible)')
            end
            obj.cur_state = obj.LatticeValue{1};
            obj.t = 1;
        end
        
        function state_ok = checkState(obj, t, state)
            % CHECKSTATE Check that state is valid for a given time
            %
            % rand_proc_object.checkState(t, state)
            %       Raise 'RandProc:InvalidState' error if t is not valid in time t
            % state_ok = rand_proc_object.checkState(t, state)
            %       No error, simply return true/false if state is
            %       valid/not
            state_ok = not(isempty(state)) && ismembertol(state, obj.LatticeUncondProb{t}, obj.tol);
            
            if nargout == 0 && not(state_ok)
                error('RandProcess:InvalidState', 'State %g is not valid at time %d', state, t)
            end
        end


    end

end
