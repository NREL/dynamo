function [ post_state_list ] = HEMS_DecisionApply( params, t, pre_state_list, decision_list )
%function [ post_state_list ] = HEMS_DecisionApply( params, t, pre_state_list, decision_list )
% Generate the post-decision state list.
% Adapted from MultiInv by BryanP
% Author: Li Wang


% Handle "zero" decision case for adpTD1
if isequal(decision_list, 0)
    post_state_list = pre_state_list;
    return
end

% Find the corresponding next states
post_state_list(:,1)= round(pre_state_list(1) + decision_list(:,1),1);%HVAC unit space 0.1

end

