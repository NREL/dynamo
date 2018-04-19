function [ next_pre_state_list ] = HEMS_RandomApply( params, t, post_state_list, uncertain_demand_list )
%function [ next_pre_state_list ] = HEMS_RandomApply( params, t, post_state_list, uncertain_demand_list )
%Adapted from MultiInv by BryanP
%Author: Li Wang

coef_1 = [0.9];%%coef of outside temp, coef of water usage, coef of EV
coef_2 = [0.1];
% pre_state_list = post_state_list;


%calculate the next pre state
next_pre_state_list(:,1) = round(coef_1(1).*post_state_list(:,1) + coef_2(1).*uncertain_demand_list(:,1),1);%HVAC unit space 0.1
%below code is to limit the state inside of the appliance range
cap_mat = cell2mat(params.appliance_range);
l_min = cap_mat(:,1)';
l_max = cap_mat(:,2)';
next_pre_state_list = max(next_pre_state_list,l_min);
next_pre_state_list = min(next_pre_state_list,l_max);

end

