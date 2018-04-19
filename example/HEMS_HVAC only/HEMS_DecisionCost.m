function [ c_list ] = HEMS_DecisionCost( params, t, pre_state_list, decision_list )
%function [ c_list ] = HEMS_DecisionCost( params, t, pre_state_list, post_state_list )
% Calculate decision cost(electricity cost).
% Adapted from MultiInv by BryanP
% Author: Li Wang
coef_ele = params.electricity_price(t);    % Price of electricity, $/KWh
 
c_list = abs(decision_list)* params.unit_cost' * coef_ele; % params.unit_cost is degree KWh/Degree F for HVAC and WH, or KWh/ 1 percent
%SOC for EV
end

