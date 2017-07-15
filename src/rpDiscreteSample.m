classdef rpDiscreteSample < RandProcess
%rpDiscreteSample Sample from diff process at each state
%
% Random process that samples from a different discrete distribution
% for each time period independant of previous time period's state
%
% Constructor: obj = rpDiscreteSample(v_list, p_list)
%
% Required inputs/properties:
%  v_list: cell row vector with one cell per time period containing a
%           numeric col vector of possible values
% Optional:
%  p_list: matching cell vector with corresponing probabilites for each
%           value. The sum of each col must equal 1. If p_list is not
%           provided, a uniform distribution across all values in the
%           corresponding v_list is assumed
%
% Additional properties (set using obj.PROPERTY):
%  None
%
% Notes:
%  - For any t greater than the max set in the ValueList, the final value
%    from the list is assumed to be held constant
%  - For non-integer times, the output is assumed to remain at the value
%    corresponding to the previous integer time step (Zero-order hold)
%
% Examples:
% % 1-D state with constant (in time) values and specified probability distribution 
% >> rng(0); s = rpDiscreteSample({[0, 1, 2, 3]'}, {[0.1, 0.5, 0.35, 0.05]'});
% >> s.sample()
% ans = 
%        2
%
%
% % 2-D state with different first period values and automatic uniform distribution
% >> s = rpDiscreteSample({[0 1; 1 2; 3 4; 5 6], [10 11; 20 21; 30 31]});
% >> [t1_vals, t1_prob] = s.dlistnext()
% 
% t1_vals =
% 
%     10    11
%     20    21
%     30    31
% 
% 
% t1_prob =
% 
%     0.3333
%     0.3333
%     0.3333
%
% >> s.cur_state()
%
% ans =
% 
%      0     1
% 
% >> s.t
% 
% ans =
% 
%      1
% 
% >> s.step()
% 
% ans =
% 
%     30    31
% 
% >> s.t
% 
% ans =
% 
%      2
% 
% >> s.step()
% 
% ans =
% 
%     20    21
% 
% >> s.step()
% 
% ans =
% 
%     10    11
% 
% >> s.t
% 
% ans =
% 
%      4
%
%
%--Range tests
% >> s = rpDiscreteSample({[0 1; 1 2; 3 4; 5 6], [10 11; 20 21; 30 31]});
% >> s.range()
%
% ans =
%
%      0     1
%      5     6
%
% >> s.range(2)
%
% ans =
%
%     10    11
%     30    31
%
% >> s.range('all')
%
% ans =
%
%      0     1
%     30    31
%
%
% see also RandProc, rpDeterministic, rpMarkov, rpLattice
%
% originally by Bryan Palmintier 2011

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  15  2017-07-14 21:35  BryanP      Remove sim and dsim, distinquish N_uniqueT from Tmax 
%  14  2017-07-14 11:10  BryanP      Refactor generally usable code to RandProc 
%  13  2017-04-25 07:42  BryanP      Allow empty state since we are time independant
%  12  2017-04-05 23:42  BryanP      Bugfixes & add checkState()
%  11  2017-04-05 22:02  BryanP      Overhaul complete:
%                                       -- t now 1 indexed
%                                       -- Removed any numbered states. Now only value-based states 
%                                       -- Removed extraneous functions including backward listing 
%                                       -- Put t as first parameter for most methods, that use t, except sample 
%  10  2016-04-05 14:53  BryanP      First step to value-only states: fixed dlistnext 
%   9  2016-xx           BryanP      Corrected time indexing (reported to dynamo 3/2017) 
%   8  2016-05-01 03:25  BryanP      Expose sample() and add support for vectorized samples
%   7  2016-04-29 00:25  BryanP      Use uniform probability if P_List not provided
%   6  2012-04-17 08:55  BryanP      Reworked cdf for sim/sample
%   5  2011-11-15 18:30  BryanP      Corrected bug with sample return order
%   4  2011-06-09 16:00  BryanP      match input shape for dnum2val and dval2num
%   3  2011-05-03 01:04  BryanP      Fixed:
%                                      - time/index in dlistnext & prev
%                                      - random sample on reset
%                                      - dlist state order for non-ascending values
%   2  2011-04-12 09:45  BryanP      Initial coding complete
%   1  2011-04-08 14:30  BryanP      Adapted from rpLattice v6


    % Internal properties
    properties (Access='protected')
        Tol = 1e-6;     % Tolerance for checking probabilities sum to 1
    end

    methods
        %% ===== Constructor & related
        function obj = rpDiscreteSample(v_list, p_list)
            % Support zero parameter calls to constructor for special
            % MATLAB situations (see help files)
            if nargin == 0
                v_list = obj.Values;
                p_list = obj.UncondProbs;
            end
            
            %-- Extract dimensional data
            %Note: N_uniqueT is the max t with a unique value list. passing
            %a t>N_uniqueT will simply return the N_uniqueT values (zero order
            %hold)
            obj.N_dim = size(v_list{1},2);

            % If probabilites not provided, assume uniform across all
            % values
            if nargin < 2
                obj.UncondProbs = cell(size(v_list));
                for t_idx = 1:length(v_list)
                    num_states = size(v_list{t_idx},1);
                    obj.UncondProbs{t_idx} = ones(num_states,1) ./ num_states;
                end

            else
                obj.UncondProbs = p_list;
            end
            
            % Store state value list
            obj.Values = v_list;


            % Loop through time for setting up additional parameters
            for t_idx = 1:obj.N_uniqueT
                %Create Cumulative distribution function for easier
                %mapping of random samples
                obj.UncondCdfs{t_idx} = cumsum(obj.UncondProbs{t_idx});

                %Check for valid probability vectors 
                %  Probability must sum to 1
                if abs(obj.UncondCdfs{t_idx}(end) - 1) > obj.Tol
                    error('RandProc:ProbSum1', 'Probabilities must sum to 1.0 for t=%d', t_idx)
                end

                %  Must be same size as value list
                if size(obj.UncondProbs{t_idx}) ~= size(obj.Values{t_idx})
                    error('rpDiscreteSample:ValProbMismatch', 'Value & Probability vectors must have equal size for time=%d', t_idx)
                end
            end
            obj.reset();
        end
        
        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.

        function [next_state_list, probability_list] = dlistnext (obj, t, cur_state)
        % DLISTNEXT List next discrete states & probabilities
        %
        % List possible next states (by number) along with conditional
        % probability P(s_{t+1} | s_t) Since the samples are independant
        % of each other (ie not path dependant), this simplifies to
        % P(s_{t+1})
        %
        % If t and/or state are not defined, the current simulation time
        % and state are assumed
            if nargin == 1
                [next_state_list, probability_list] = obj.dlistnext(obj.t, obj.cur_state);
                return
            end
            
            % adapt time to give valid index
            t = min(floor(t), obj.N_uniqueT);
            
            % Allow empty state since we are time independant
            if not(isempty(cur_state))
                obj.checkState(t, cur_state);
            end

            %if we get here, we know the state is valid
            %
            % For the discrete sample, each sample is independant, so
            % just return the next periods value and sample options
            %
            %Note: t+1 is right b/c want next state
            [next_state_list, probability_list] = obj.state_info(t+1);
        end

    end
 
    methods (Access = protected)
        function state_list = conditionalSample(obj, N, t, cur_state) %#ok<INUSD>
        %CONDITIONALSAMPLE draw state samples for the specified state
        %
        % Since we don't have any time correlation, this is just a simple
        % wrapper
            state_list = obj.sample(N, t);
        end
    end
    
end
