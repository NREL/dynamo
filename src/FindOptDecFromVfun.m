function [decision, dec_contrib, post_state, forward_val] = ...
                        FindOptDecFromVfun(problem, t, pre_state, vfun, vfun_approx_params)
% FINDOPTDECFROMVFUN Find optimal decision using value function
%
% [decision, dec_contrib, post_state, forward_val] = ...
%                        FindOptDecFromVfun(problem, t, pre_state, vfun)
% [decision, dec_contrib, post_state, forward_val] = ...
%                        FindOptDecFromVfun(problem, t, pre_state, vfun, vfun_approx_params)

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   7  2017-04-02 07:03  BryanP      Updated for standard problem definition 
%   6  2016-04-15 02:03  BryanP      Error checks for missing fields
%   5  2012-07-06 16:55  BryanP      Added time to PostToVfun calls
%   4  2012-05-08 14:05  BryanP      Renamed from PreToPost to ApplyDscn
%   3  2012-04-21 13:20  BryanP      Rely on faLocalRegr for auto-expand of neighborhood
%   2  2012-04-20 06:45  BryanP      Separate postdec from vfun state spaces
%   1  2012-04-16 19:15  BryanP      Extracted from parTD1 v1

    if nargin < 5 || isempty(vfun_approx_params)
        vfun_approx_params = {};
    end

    %Extract list of all possible decisions
    decision_set = problem.fDecisionSet(problem.params, pre_state, t);
    decision_list = decision_set.as_array();
    
    %Compute corresponding costs
    possible_dec_contribs = problem.fDecisionCost(problem.params, t, pre_state, decision_list);
    %And resulting post-decision states
    post_state_list = problem.fDecisionApply(problem.params, t, pre_state, decision_list);
    

    if not(isfield(problem, 'fMapState2Vfun')) || isempty(problem.fMapState2Vfun)
        vfun_state_list = post_state_list;
    else
        vfun_state_list = problem.fMapState2Vfun(problem, post_state_list, t);
    end

    %Find post-decision value for all possible decisions
    future_val_list = approx(vfun, vfun_state_list, vfun_approx_params{:});

    total_cost_list = possible_dec_contribs + future_val_list;

    %Find optimal decision
    [~, best_idx] = max(total_cost_list);

    decision = decision_list(best_idx, :);
    dec_contrib = possible_dec_contribs(best_idx);
    post_state = post_state_list(best_idx, :);
    forward_val = future_val_list(best_idx);
end
