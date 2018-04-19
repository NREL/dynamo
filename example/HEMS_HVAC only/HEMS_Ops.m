function [ ops_contribs ] = HEMS_Ops( params, t, next_pre_state_list, this_decision, random_outcome_list )
%function [ after_rand_ops ] = HEMS_Ops( params, t, next_pre_state_list, this_decision, random_outcome_list )
% Calculate the discomfort for each applicance. NOTE: There is No discomfort
% for EV since it is only calculated in final state(Terminal Value)
% Adapted from MultiInv by BryanP
% Author: Li Wang

%coef_cost = [-0.06;-0.019;-0];% Dissatisfation coefficient for deviation from desired room temp, water temp and EV SOC
coef_cost = [-0.477];
ops_contribs = (abs(next_pre_state_list - params.desired_value))*coef_cost;


end

