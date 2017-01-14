function [results] = dpBI(problem, varargin)
% dpBI Backward Induction for DP and intialization
%
% [results] = dpBI(problem, verbose, true, discount_rate, 0.9)
% [results] = dpBI(problem, dp_opt, struct('verbose', true))
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
%   2  2016-10-27          dkrishna   Add core algorithm
%   1  2016-03-07          dkrishna   Pseudo-code


p = inputParser;
p.addRequired('problem');
p.addParameter('verbose', false);
p.addParameter('discount_rate', false);
p.addOptional('dp_opt', struct());
p.parse(problem, varargin{:});

problem = p.Results.problem;

verifyProblemStruct(problem);

dp_opt = parse_dp_opt(p.Results.dp_opt, p);

if dp_opt.verbose
    fprintf('Backward Induction DP\n')
end

cell_values = cell(1, problem.n_periods);
cell_policy = cell(1, problem.n_periods);

for period = problem.n_periods:-1:1
    states = problem.state_set{period}.as_array();
    n_states = length(states);

    values = -Inf * ones(n_states, 1);
    policy = NaN * ones(n_states, 1); %action (orders) = f(s,t), no action in final state
    policy = num2cell(policy);

    if period == problem.n_periods

        [v] = problem.fTerminalValue(problem.params, period, states);
        values(1:n_states) = v;
        
    else
        for s_idx = size(states, 1)
            
            s = states(s_idx, :);

            decision_set = problem.fDecisionSet(problem.params, s, period);
            [decisions] = decision_set.as_array();

            if ~ isempty(decisions)

                n_decisions = length(decisions);

                contribution = zeros(n_decisions, 1);
                for d = decisions(:)'
                    % post_decision_state = problem.fDecisionApply(problem.params, s, d, period);
                    value_from_decision = problem.fDecisionCost(problem.params, s, d, period);
                    
                    n_random_process = length([problem.random_items{:}]);
                    for rp = [problem.random_items{:}]
                        
                        [value_list, state_n_list, prob] = rp.dlistnext(s, period);

                        random_process_value = prob*problem.fRandomProcessValue(value_list);
                        value_from_decision = value_from_decision + random_process_value;
                    end

                    % TODO combine to a a realizations and probablity array
                    % TODO different value function for every time period
                    % Integer range from real
                    % real range from integer
                    contribution(d+1) = value_from_decision;
                    % [uncertainty_value, problem.params] = problem.fUncertaintyApply(s, d, period, cell_values{period+1}, problem.params);
                end

                % TODO See if FindOptimalDecision is available

                [best_contribution, best_contribution_index] = max(contribution);

                if best_contribution > values(s)
                    values(s) = best_contribution;
                    policy{s} = best_contribution_index-1;
                end
            end
        end
    end

    cell_values{period} = values;
    cell_policy{period} = policy;
end


results.policy = cell_policy;
% results.firstPeriodDecision = firstPeriodDecision;
% results.firstPeriodObjectiveFunction = firstPeriodObjectiveFunction;
results.opts = dp_opt;
results.values = cell_values;

end % Main Function

function dp_opt = parse_dp_opt(dp_opt, p)
% parse_dp_opt Handle dynamic programming option parsing
%
% dp_opt = parse_dp_opt(struct('verbose', true))
%
% Required
%   dp_opt
%

if ~ isfield(dp_opt, 'verbose')
   dp_opt.verbose = p.Results.verbose;
end
if ~ isfield(dp_opt, 'discount_rate')
   dp_opt.discount_rate = p.Results.discount_rate;
end


end

