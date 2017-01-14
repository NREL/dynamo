function [start_decision, start_post] = MultiInvDecisionSample(params, t, pre_s_list, ~)
% Hacked together function to support reasonable decision samples for
% adpTD1. Provides one sample and resulting post decision state per each
% pre_s_list item 

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2016-05-01 04:11  BryanP      Initial Code 
num_s = size(pre_s_list,1);
start_decision = nan(num_s, params.dims.decision);
start_post = nan(num_s, params.dims.post_state);

for s_idx = 1 : size(pre_s_list,1)
    % ID number of possible decisions
    num_desc_for_s = MultiInvDecision(pre_s_list(s_idx, :), 0, t, [], params);
    % Sample across this number and extract the decision of interest
    start_decision(s_idx, :) =  MultiInvDecision(pre_s_list(s_idx, :), randi(num_desc_for_s), t, [], params);
    % And find corresponding post decision state
    start_post(s_idx, :) = MultiInvApplyDecision(params, pre_s_list(s_idx, :), start_decision(s_idx, :), t);
end
    
    