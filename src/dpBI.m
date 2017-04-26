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
%   6  2017-04-26 17:32    BryanP     BUGFIX: correct operations costs (now Working!) 
%   5  2017-04-26 05:32    BryanP     Working? 
%   4  2017-04-25 07:00    bpalmint   WIP: mostly through uncertainty management 
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
policy = cell(1, problem.n_periods);  %Contains optimal decisions for each state


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
% Explanation: The goal here is to run th core backward induction algorithm
% of looping over all uncertainty for all decisions for all states. A
% representative basic algorithm for this would be:  
%
%>>>   list all pre_states for t
%>>>   for each pre_state
%>>>     compute before_decision operations cost
%>>>     list all possible decisions
%>>>     for each decision
%>>>       compute decision cost
%>>>       determine the resulting post decision state
%>>>       compute after decision operations costs
%>>>       list all uncertainty outcomes and corresponding probabilities
%>>>       for each uncertainty outcome
%>>>         compute uncertainty contribution
%>>>         determine the resulting next pre-decision state
%>>>         compute after uncertainty operations costs
%>>>         extract value from next (t+1) pre-decision state
%>>>         compute post_random_val as sum( uncertainty_contribution, after_random_ops, (1-disc_rate) * next_pre_value) 
%>>>       compute E[random_value] as sum( post_random_val * probability )
%>>>       compute total_decision_value as sum(after_dec_ops_cost, expected_random_value)
%>>>     find optimal_decision as argmin(total_decision_values) 
%>>>     store optimal_decision in policy array for pre_state and t
%>>>     compute total_expected_value for pre-decision statate as sum(total_decision_value, before_decision_ops) 
%>>>     store total_expected_value in value array for pre_state and t
%
% However...

for t = problem.n_periods:-1:1
    
    if dp_opt.verbose
        fprintf('    T=%d:', t)
    end
    
    %>>>   list all pre_states for t
    %extract (pre-decision) states
    pre_state_list{t} = problem.state_set{t}.as_array();
    n_pre_states = size(pre_state_list{t}, 1);

    %>>>   for each pre_state
    %>>>     compute before_decision operations cost (vectorized)
    if not(isempty(problem.fOpsBeforeDecision))
        % Compute unique ops costs, using internal Ops loop
        before_dec_ops = problem.fOpsBeforeDecision(params_only, t, pre_state_list{t});
    else
        before_dec_ops = zeros(n_pre_states, 1);
    end

    %Initialize optimal policy (decision) storage
    first_decision_set = problem.fDecisionSet(problem.params, t, pre_state_list{t}(1, :));
    n_decision_dims = first_decision_set.N_dim;
    policy{t} = NaN(n_pre_states, n_decision_dims);
    values{t} = NaN(n_pre_states, 1);
        
    %Find all possible decisions for all possible pre-states and produce
    %set of resulting post_states
    for s = 1:n_pre_states
        
        DisplayProgress(dp_opt.verbose, s)
        
        this_pre_state = pre_state_list{t}(s, :);
        
        %>>>     list all possible decisions
        decision_set = problem.fDecisionSet(problem.params, t, this_pre_state);
        decision_list = decision_set.as_array();
        n_decisions = size(decision_list, 1);
                
        %>>>     for each decision
        %>>>       compute decision cost (vectorized)
        decision_cost_list = problem.fDecisionCost(problem.params, t, this_pre_state, decision_list);
        %>>>       determine the resulting post decision state (vectorized)
        post_state_list = problem.fDecisionApply(problem.params, t, this_pre_state, decision_list);
        
        %>>>       compute after decision operations costs (vectorized)
        if not(isempty(problem.fOpsAfterDecision))
            % Compute unique ops costs, using internal Ops loop
            after_dec_ops = problem.fOpsAfterDecision(problem.params, t, post_state_list, decision_list);
        else
            after_dec_ops = 0;
        end
    
        %actually start decision loop for non-vecorizable pieces
        expected_random_val = NaN(n_decisions, 1);
        for d = 1:n_decisions
            this_decision = decision_list(d, :);
            this_post_state = post_state_list(d, :);
            
            %>>>       list all uncertainty outcomes and corresponding probabilities
            % extract random portion of the state:
            if not(isempty(problem.random_state_map))
                cur_rand_state = utilRandStatefromState(this_post_state, problem.random_state_map);
            else
                cur_rand_state = {};
            end
            [random_outcome_list, probability_list] = RandSetNextJoint(problem.random_items, t, cur_rand_state);
            
            %>>>       for each uncertainty outcome
            %            NOTE: fully vectorized so no loop
            %>>>         compute uncertainty contribution
            uncertainty_contrib = problem.fRandomCost(problem.params, t, this_post_state, random_outcome_list);
            %>>>         determine the resulting next pre-decision state
            next_pre_state_list = problem.fRandomApply(problem.params, t, this_post_state, random_outcome_list);
            %>>>         compute after uncertainty operations costs
            if not(isempty(problem.fOpsAfterRandom))
                after_rand_ops = problem.fOpsAfterRandom(problem.params, t, next_pre_state_list, this_decision, random_outcome_list);
            else
                after_rand_ops = 0;
            end
            
            %>>>         extract value from next (t+1) pre-decision state
            next_pre_value = zeros(size(uncertainty_contrib));
            [valid_state_mask, next_state_map] = ismember(next_pre_state_list, pre_state_list{t+1}, 'rows');
            if not(all(valid_state_mask))
                warning('dbBI:state_not_found', 'Some post uncertainty states not found for post-state: [%s] at t=%d. Setting values to -Inf', ...
                    sprintf('%g ', this_post_state), t)
                next_pre_value(not(valid_state_mask)) = -Inf;
            end
            next_pre_value(valid_state_mask) = values{t+1}(next_state_map);
            
            %>>>         compute post_random_val as sum( uncertainty_contribution, after_random_ops, (1-disc_rate) * next_pre_value) 
            post_random_value = uncertainty_contrib + after_rand_ops + (1-problem.discount_rate) * next_pre_value;

            %>>>       compute E[random_value] as sum( post_random_val * probability )
            expected_random_val(d) = sum(post_random_value .* probability_list);
            
        end
        %>>>       compute decision_value as sum(after_dec_ops_cost, expected_random_value, decision_cost_list)
        total_decision_value = after_dec_ops + expected_random_val + decision_cost_list;
        
        %>>>     find optimal_decision as argmin(decision_values)
        optimal_total_decision_value = max(total_decision_value);
        optimal_decision_idx = find(total_decision_value == optimal_total_decision_value, 1, 'first');
        optimal_decision = decision_list(optimal_decision_idx, :);

        %>>>     store optimal_decision in policy array for pre_state and t
        policy{t}(s,:) = optimal_decision;
        
        %>>>     compute total_expected_value for pre-decision statate as sum(total_decision_value, before_decision_ops) 
        total_expected_value = optimal_total_decision_value + before_dec_ops(s);
        %>>>     store total_expected_value in value array for pre_state and t
        values{t}(s) = total_expected_value;

    end

    if dp_opt.verbose
        fprintf('Done: %d states\n', n_pre_states)
    end

end

results.dpbi_policy = policy;
% results.firstPeriodDecision = firstPeriodDecision;
% results.firstPeriodObjectiveFunction = firstPeriodObjectiveFunction;
results.dp_opts = dp_opt;
results.dpbi_values = values;
results.pre_state_list = pre_state_list;

%% ===== Clean-up =====
%Reset Auto-parallel state
if not(isempty(cache_par_auto))
    % when non-empty, we already created ps
        ps.Pool.AutoCreate = cache_par_auto;
end


end % Main Function
