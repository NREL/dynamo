function [post_state_list] = ...
    InvDecisionApply(params, t, pre_state_list, decision_list) %#ok<INUSL>
% INVDECISIONAPPLY Applies the provided decision(s) to the current (pre-decision) state 
% 
% Produces a (set of) new post decision state(s).

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-09-23        BryanP      Copied from MultiInvDecisionApply v10

% Handle "zero" decision case for adpTD1
if isequal(decision_list, 0)
    post_state_list = pre_state_list;
    return
end

% Find the corresponding next states
post_state_list = bsxfun(@plus, pre_state_list, decision_list);

