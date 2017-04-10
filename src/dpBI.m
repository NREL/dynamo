function [results] = dpBI(problem, varargin)
% dpBI Backward Induction for DP and intialization
%
% [results] = dpBI(problem, dp_opt_structure)
% [results] = dpBI(problem, 'verbose', true)
%
% Required problem attributes
%   problem
%
% Optional attributes
%   discount_rate
%   verbose
%
% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- -------------------------------------
%   3  2017-04-09 23:00    bpalmint   WIP: defaults, setup, pre-decision operations cost, parrallel decisions 
%   2  2016-10-27          dkrishna   Add core algorithm
%   1  2016-03-07          dkrishna   Pseudo-code

%% ====== Handle Inputs =====
%--- Handle DP options
%TODO: optional support for varargin list rather than 
if nargin < 2 || isempty(varargin)
    varargin = struct([]);
end

dp_defaults = {
                %general DP options used here
                'verbose'               50
                'parallel'              false
                'fix_rand'              false
                'fix_rand_is_done'      false
               };

dp_opt = DefaultOpts(varargin, dp_defaults);

%--- Handle problem setup
verifyProblemStruct(problem);

if dp_opt.verbose
    fprintf('Backward Induction DP\n')
end

%initialize output storage
% TODO: consider directly storing in results structure
values = cell(1, problem.n_periods + 1);
pre_state_list = cell(1, problem.n_periods + 1);
policy = cell(1, problem.n_periods);


%% ====== Standardized Setup =====
% Configure random numbers, including setup consistent stream if desired
dp_opt = utilRandSetup(dp_opt);

% Initialize parallel pool if required (and suppress parallel pool
% initializtion if parallel off)
[cache_par_auto, ps] = utilParSetup(dp_opt);    

% Add additional required functions as needed
problem.fOptimalDecision = utilFunctForProblem(problem, 'fOptimalDecision', @FindOptDecFromVfun);
problem.fRandomSample = utilFunctForProblem(problem, 'fRandomJoint', @RandSetNextJoint);

%% ====== DP Backward Induction Algorithm =====
%% --- terminal period (T+1) ---
t = problem.n_periods + 1;
% TODO: Use problem specific sample function if one is defined
if dp_opt.verbose
    fprintf('    T=%d (terminal period): ', t)
end

pre_state_list{t} = problem.state_set{t}.as_array();

% -- Compute contribution for sampled states
% Note: Memoized operations in terminal period must be handled by user. Typically
% by calling a memoized ops function in fTerminalValue and storing the Ops
% table in a handle derived class.
values{t} = problem.fTerminalValue(problem.params, t, pre_state_list{t});

if dp_opt.verbose
    fprintf('Done\n')
end

%% --- earlier periods ---
for t = problem.n_periods:-1:1
    
    if dp_opt.verbose
        fprintf('    T=%d:', t)
    end
    
    %Cache data needed for parfor
    fn_decision_set = problem.fDecisionSet;
    fn_decision_apply = problem.fDecisionApply;
    params_only = problem.params;

    %extract (pre-decision) states
    pre_state_list{t} = problem.state_set{t}.as_array();
    
    %Setup decision-related storage
    n_pre_states = size(pre_state_list{t}, 1);
    post_state_list = cell(n_pre_states, 1);
    decision_contrib = cell(n_pre_states, 1);
    cost_post_dec_ops = cell(n_pre_states, 1);

    % compute before decision operations costs
    if not(isempty(problem.fOpsBeforeDecision))
        % Compute unique ops costs, using internal Ops loop
        before_dec_ops = problem.fOpsBeforeDecision(params_only, t, pre_state_list{t});
    else
        before_dec_ops = zeros(n_pre_states, 1);
    end

    parfor s = 1:n_pre_states
        %Extract possible valid decisions
        decision_set = fn_decision_set(params_only, t, pre_state_list{t}(s, :)); %#ok<PFBNS>
        decision_list = decision_set.as_array();
        
        %Store corresponding resulting states
        post_state_list{s} = fn_decision_apply(params_only, t, pre_state_list{t}(s, :), decision_list); %#ok<PFBNS>
    end

%>>>>>>>>>>>>>>>> BP EDIT MARKER
    
    n_post_states = size(post_state_list, 1);
    uncertainty_list = cell(n_post_states, 1);

    % Note: since we typically sample multiple times per post decision
    % state, each result is stored as a cell array
    % Cache only required structure pieces for use in the following set of
    % parfor loops
    % Note: some values, e.g. params already cached
    fn_random_apply = problem.fRandomApply;
    fn_random_cost = problem.fRandomCost;
    fn_random_sample = problem.fRandomSample;
    fn_optimal_decision = problem.fOptimalDecision;
    %Note: Decision costs not included b/c assume internal loop for computing operations cost (with memoization?) 
    n_periods = problem.n_periods;
    rand_per_post_state = 10;
    
    for post_idx = 1:n_post_states
        if dp_opt.verbose
            DisplayProgress(dp_opt.verbose, post_idx)
        end
        
        %>>>     sample uncertainty and store change
        % Sample Random outcomes to get to next pre-states
                
        rand_list = problem.random_items{t};
        if length(rand_list) > 1
            warning('dpBI:NotImplemented', 'Currently Random Items can only be lenght 1 at each time') 
        end
        
        rand_list.dlistnext(post_state_list(post_idx,:), t)
        
        display('Stop')
        
        next_pre_list{post_idx} = fn_random_apply(params_only, t, post_state_list(post_idx,:), rand_list{rand_idx}.Vlist);
        
    end
    
%     for s = 1:length(state_list)
%         %Extract valid decisions
%         decision_set = fn_decision_set(params_only, t, states(s, :)); %#ok<PFBNS>
%         %And sample these choices
% 
%         
%         ds = decision_set.as_array();
%         decision_set
%         decision_list{s} = ds{s};
%         %Build a piece of the post_decision sample
%         state_list{s} = fn_decision_apply(params_only, t, state_list(s, :), decision_list{s}); %#ok<PFBNS>
%     end
%     
%     % Finally reconstruct the possibly parallel pieces into a full list
%     state_list = cell2mat(state_list);
%     decision_list = cell2mat(decision_list);
    

end


% results.policy = cell_policy;
% % results.firstPeriodDecision = firstPeriodDecision;
% % results.firstPeriodObjectiveFunction = firstPeriodObjectiveFunction;
% results.opts = dp_opt;
% results.values = cell_values;

%% ===== Clean-up =====
%Reset Auto-parallel state
if not(isempty(cache_par_auto))
    % when non-empty, we already created ps
        ps.Pool.AutoCreate = cache_par_auto;
end


end % Main Function



% 
%     values = -Inf * ones(n_states, 1);
%     policy = NaN * ones(n_states, 1); %action (orders) = f(s,t), no action in final state
%     policy = num2cell(policy);
% 
%     if t == problem.n_periods
% 
%         [v] = problem.fTerminalValue(problem.params, t, states);
%         values(1:n_states) = v;
%         
%     else
%         for s_idx = size(states, 1)
%             
%             s = states(s_idx, :);
% 
%             decision_set = problem.fDecisionSet(problem.params, s, t);
%             [decisions] = decision_set.as_array();
% 
%             if ~ isempty(decisions)
% 
%                 n_decisions = length(decisions);
% 
%                 contribution = zeros(n_decisions, 1);
%                 for d = decisions(:)'
%                     % post_decision_state = problem.fDecisionApply(problem.params, s, d, period);
%                     value_from_decision = problem.fDecisionCost(problem.params, s, d, t);
%                     
%                     n_random_process = length([problem.random_items{:}]);
%                     for rp = [problem.random_items{:}]
%                         
%                         [value_list, state_n_list, prob] = rp.dlistnext(s, t);
% 
%                         random_process_value = prob*problem.fRandomProcessValue(value_list);
%                         value_from_decision = value_from_decision + random_process_value;
%                     end
% 
%                     % TODO combine to a a realizations and probablity array
%                     % TODO different value function for every time period
%                     % Integer range from real
%                     % real range from integer
%                     contribution(d+1) = value_from_decision;
%                     % [uncertainty_value, problem.params] = problem.fUncertaintyApply(s, d, period, cell_values{period+1}, problem.params);
%                 end
% 
%                 % TODO See if FindOptimalDecision is available
% 
%                 [best_contribution, best_contribution_index] = max(contribution);
% 
%                 if best_contribution > values(s)
%                     values(s) = best_contribution;
%                     policy{s} = best_contribution_index-1;
%                 end
%             end
%         end
%     end
% 
%     cell_values{t} = values;
%     cell_policy{t} = policy;
