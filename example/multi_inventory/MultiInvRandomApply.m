function [next_pre_state_list] = ...
    MultiInvRandomApply(params, t, post_state_list, uncertain_demand_list) %#ok<INUSL>
% MULTIINVRANDOMAPPLY Applies the provided uncertainty(ies) to the current (post-decision) state 
% 
% Produces a (set of) new pre decision state(s) to be used during the next time period. 

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2016-10-27 13:31  BryanP      Initial Code adapted from MultiInvDecisionApply v7 

% Find the corresponding next states by subtracting off demand
next_pre_state_list = bsxfun(@minus, post_state_list, uncertain_demand_list);
% And prevent any negative inventory states
next_pre_state_list = max(0, next_pre_state_list);
