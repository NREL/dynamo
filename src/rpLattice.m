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
%  CondProb:   vector of corresponding coefficients
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
% Examples: (Note: extensive, comprehensive debugging in testrpLattice unit test)
% >> lattice_object = rpLattice(5, [0.5  1.2 1.5 ]', [0.15 0.5 0.35]', 3)
% 
% lattice_object = 
% 
%   rpLattice with properties:
% 
%            Start: 5
%             Coef: [3×1 double]
%         CondProb: [3×1 double]
%              Tol: 1.0000e-04
%                t: 1
%        cur_state: 5
%             name: ''
%     pt_dim_names: {}
%       SampleType: 'rand'
%            N_dim: 1
%     DiscreteMask: []
% 
%
% see also rpTransMatrix, RandProc, rpDiscreteSample, rpMarkov, rpBasic
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  17  2017-07-16 00:13  BryanP      Use standardized conditionatlSample() in RandProcess 
%  16  2017-07-15 22:13  BryanP      Use RandProcess reset() b/c no longer supports multiple starting conditions 
%  15  2017-07-15 11:27  BryanP      BUGFIX: only build lattice to Tmax (artifact of t_start = 1) 
%  14  2017-07-15 weehr  BryanP      Overhaul: t_start = 1, RandProcess streamlining, debugging with testrpLattice 
%  13  2017-07-14 05:48  BryanP      Only support value states 
%  12  2017-07-13 21:17  BryanP      Initial Update for new RandProc format
%  11  2012-04-17 08:55  BryanP      Reworked cdf for sim/sample
%  10  2012-01-25 13:45  BryanP      Allow blank constructor input for loading from a file
%   9  2012-01-11 17:00  BryanP      Added constant non-walk for t>LatticeTmax
%   8  2011-11-09 23:40  BryanP      Added tolerance to constructor
%   7  2011-06-09 16:00  BryanP      match input shape for dnum2val and dval2num
%   6  2011-03-12 12:45  BryanP      Reordered set.Coef & CondProb to work with save/load
%   5  2011-01-04 22:00  BryanP      Complete, working version
%   4  2010-12-23 15:45  BryanP      Implemented buildLattice & dlist*
%   3  2010-12-21 21:20  BryanP      Made LatticeTmax required in constructor, etc
%   2  2010-12-15 23:30  BryanP      Added t_max for dlist(... 'all')
%   1  2010-12-16 14:30  BryanP      Adapted from rpList v6

    properties
        Start = [];   % initial value
        Coef = [];    % vector of coefficients
        CondProb = [];    % vector of probabilities
    end

    % Internal properties
    properties (Access='protected')
        ConditionalCdf = [];          %precomputed conditional probability (based on CondProb)
        LatticeTmax = 1;
        
        %Flags
        SkipLatticeBuild = false    %Wait till all values set before building the lattice
        SuppressLengthError = false %Prevent coef & prob length checks when changing both
    end

    methods (Access = protected)
        %Reset the stored Lattice
        function clearStoredLattice(obj)
            obj.Values = {obj.Start};
            obj.UncondProbs = {1};
            obj.UncondCdfs = {1};
        end

        %Build & cache the lattice
        function buildLattice(obj)
            if obj.SkipLatticeBuild
                return
            end
            % -- Build the value lattice.
            %compute conditional cdf here to avoid parameter
            %cross-reference issues
            obj.ConditionalCdf = cumsum(obj.CondProb);

            %Starting with the initial value
            obj.Values = {obj.Start};
            obj.UncondProbs = {1};
            obj.UncondCdfs = {1};
            %Then cycling through future time scenarios
            for idx=2:(obj.LatticeTmax)
                % Compute next lattice state by
                % 1) multiplying the previous set of states by the
                % coefficient matrix to form an array & compute
                % corresponding unconditional probability (for sampling)
                %Implementation Note: have to put coef first, otherwise,
                %the first multiplication creates a row vector and unique
                %does not change it back to a column vector.
                obj.Values{idx} = obj.Coef * obj.Values{idx-1}';
                obj.UncondProbs{idx} = obj.CondProb * obj.UncondProbs{idx-1}';

                % 2) rounding these values to make sure that states merge
                % Note: the RoundTo function is part of the adp toolbox
                obj.Values{idx} = RoundTo(obj.Values{idx}, obj.Tol);

                % 3) converting the resulting rectangular matrices to
                % column vectors.
                %
                % Note: in reshape, the [] will compute the required number
                % of rows.
                obj.Values{idx} = reshape(obj.Values{idx},[],1);
                obj.UncondProbs{idx} = reshape(obj.UncondProbs{idx},[],1);

                % 4) Using unique() to return an ordered vector of unique
                % values.
                [obj.Values{idx}, junk, idx_list]  = unique(obj.Values{idx}); %#ok<ASGLU>

                % 5) Use the indices (idx_list) to locate the corresponding
                % probabilities and then sum the combined probabilities
                % using accumarray
                obj.UncondProbs{idx} = accumarray(idx_list, obj.UncondProbs{idx});

                % 6) Normalize the probability so we are sure to sum to 1
                obj.UncondProbs{idx} = obj.UncondProbs{idx}/sum(obj.UncondProbs{idx});
                
                % 7) Compute corresponding unconditional cdf
                obj.UncondCdfs{idx} = cumsum(obj.UncondProbs{idx});
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
            if nargin == 0
                return
            end
            
            if nargin < 5 || isempty(tol)
                tol = 0.0001;
            end
            
            if nargin < 4 || isempty(t_max)
                t_max = 2;
            end
            
            obj.setparams(start, coef, prob, t_max, tol);
            obj.reset();
        end

        function setparams(obj, start, coef, prob, t_max, tol)
        %SETPARAMS allow easy changing of all lattice parameters at once
            %Prevent building the stored lattice until required values set
            obj.SkipLatticeBuild = true;
            %Prevent errors from set functions for length mis-matches
            obj.SuppressLengthError = true;

            if nargin > 4
                obj.LatticeTmax = t_max;
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
            obj.CondProb = prob;

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
            % Note: we do this after assigning obj.CondProb to enable
            % successful saves & loads.
            if not(obj.SuppressLengthError) %#ok<MCSUP>
                if length(c) ~= length(obj.CondProb) %#ok<MCSUP>
                        %Allow empty for loading Lattice objects from *.mat files
                        if not(isempty(c)) && not(isempty(obj.CondProb)) %#ok<MCSUP> OK b/c we explicitly check for empty (unset) values
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
        function set.CondProb(obj, p)
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
            obj.CondProb = p;

            % Ensure prob length matches coef length, unless this check is
            % suppressed b/c setting both
            %
            % Note: we do this after assigning obj.CondProb to enable
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
        
        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.

        function [state_list, prob] = dlistprev (obj, state_in, t )
        % DLISTPREV List previous discrete states & probabilities
        %
        % List possible previous states along with conditional
        % probability P(s_{t-1} | s_t)
        %
        % If t is not defined, the current simulation time is assumed
        %
        % Note: provided for user convenience 
            if nargin < 3
                if nargin < 2
                    state_in = obj.cur_state;
                end
                [state_list, prob] = obj.dlistprev(state_in, obj.t);
                return
            elseif t <= 1
                error('RandProcess:InvalidTime', 'Only t>1 valid for rpLattice')
            end

            %find a valid time for state lookup, by limiting to LatticeTmax
            t_lookup = min(t, obj.LatticeTmax);

            if isempty(state_in) ...
                    || not(all(ismembertol(state_in,obj.Values{t_lookup}, obj.Tol)))
                error('RandProcess:InvalidState', 'State %g not valid at time=%d', state_in, t)
            else
                %if we get here, we know the state is valid for this time

                if t > obj.LatticeTmax
                    state_list = state_in;
                else

                    % Build previous value list by reversing the multiplication
                    % by coef required to get there
                    state_list = RoundTo(state_in ./ obj.Coef, obj.Tol);

                    % For values near the "edge" not all of the possible priors
                    % from this division are actually valid states during the
                    % previous period, so limit our seach to those that are.
                    % Note: that the index t, corresponds to t-1 b/c 1-indexed
                    [state_list, coef_used, states] = intersect(state_list, obj.Values{t-1});
                end

                if nargout > 1
                    if t > obj.LatticeTmax
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
                        prob = obj.CondProb(coef_used) .* obj.UncondProbs{t}(states);
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
            if nargin < 2 || isempty(state_in)
                state_in = obj.cur_state;
            end
            
            if nargin < 3 || isempty(t)
                t = obj.t;
            end
            
            obj.checkState(t, state_in);
            %if we get here, we know the state is valid for this time

            if t >= obj.LatticeTmax
                state_list = state_in;
            else
                % Build next value list by multiplying
                % by coef required to get there
                state_list = RoundTo(state_in .* obj.Coef, obj.Tol);
            end

            if nargout > 1
                %Here the probability is easy... it is either
                if t >= obj.LatticeTmax
                    %one or
                    prob = 1;
                else
                    %our probability vector
                    prob = obj.CondProb;
                end
            end
        end
    end

end
