function [post_state_list] = ...
    SimpleStoragePvDecisionApply(params, t, pre_state_list, decision_list) %#ok<INUSL>
% SIMPLESTORAGEPVDECISIONAPPLY Applies the provided decision(s) to the current (pre-decision) state 
% 
% Produces a (set of) new post decision state(s).  

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-07-18 07:09  BryanP      Initial Code adapted from SimpleStoragePvRandomApply v1 

%Start with pre-decision list
post_state_list = pre_state_list;
if size(post_state_list, 1) == 1
    post_state_list = repmat(post_state_list, size(decision_list, 1), 1);
end
% Adjust storage power rating and capacity
post_state_list(:,2:end) = post_state_list(:,2:end) + decision_list;
% And prevent any negative quantities
assert(all(all(post_state_list >= 0)),'SimpStorePv:NegState', 'Unexpected Negative state value')
