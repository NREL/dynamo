function [decision, dec_contrib, post_state, forward_val] = ...
                        FindOptDecFromVfun(problem, t, pre_s, vfun, adp)
% FINDOPTDECFROMVFUN Find optimal decision using value function
%
% [decision, dec_contrib, post_state, forward_val] = ...
%                        FindOptDecFromVfun(problem, t, pre_s, vfun, adp)
%
% problem must contain function pointers (as fields)
%   fDecision:      d_list = problem.fDecision(problem, pre_s, t);
%   fApplyDscn

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2012-04-16 19:15  BryanP      Extracted from parTD1 v1
%   2  2012-04-20 06:45  BryanP      Separate postdec from vfun state spaces
%   3  2012-04-21 13:20  BryanP      Rely on faLocalRegr for auto-expand of neighborhood
%   4  2012-05-08 14:05  BryanP      Renamed from PreToPost to ApplyDscn
%   5  2012-07-06 16:55  BryanP      Added time to PostToVfun calls
%   6  2016-04-15 02:03  BryanP      Error checks for missing fields


    d_list = problem.fDecision(problem, pre_s, t);
    [post_s_list, possible_dec_contribs] = problem.fApplyDscn(problem, pre_s, d_list, t);
    if not(isfield(problem, 'fPostToVfun')) || isempty(problem.fPostToVfun)
        vfun_state_list = post_s_list;
    else
        vfun_state_list = problem.fPostToVfun(problem, post_s_list, t);
    end

    %Find post-decision value for all possible decisions (in t+1 money)
    fut_val = approx(vfun, vfun_state_list, adp.vfun_approx_params{:});

    total_cost = possible_dec_contribs + (1-problem.disc_rate) * fut_val;

    %Find optimal decision
    [~, best_idx] = max(total_cost);

    decision = d_list(best_idx, :);
    dec_contrib = possible_dec_contribs(best_idx);
    forward_val = fut_val(best_idx);
    post_state = post_s_list(best_idx, :);
end
