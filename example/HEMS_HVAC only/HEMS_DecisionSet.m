function [ decision_set, params  ] = HEMS_DecisionSet( params, t, pre_decision_state )
% function [ decision_set, params  ] = HEMS_DecisionSet( params, t, pre_decision_state )
% Generate the decision set.
% Adapted from MultiInv by BryanP
% Author: Li Wang

    max_change = [0];
    min_change = [-4]; % Here is the decision set for each applicance, the decision can take any integer value between the min_change and max_change


% Below code is to limit the decision set inside of the appliance range.
cap_mat = cell2mat(params.appliance_range);
r_min = cap_mat(:,1)';
r_max = cap_mat(:,2)';
v_max = min(pre_decision_state+ max_change,r_max);
v_min = max(pre_decision_state+ min_change,r_min);
v_max = (v_max - pre_decision_state)./params.unit_space;%number of items
v_min = (v_min - pre_decision_state)./params.unit_space;

decision_set = setCombinWithLimits('',  true, params.unit_space, params.total_space,v_max,-100,v_min);
end

