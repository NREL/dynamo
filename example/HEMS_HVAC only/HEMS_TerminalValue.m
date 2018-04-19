function [ v, params ] = HEMS_TerminalValue( params, t, s_list )
% function [ v, params ] = HEMS_TerminalValue( params, t, s_list )
% Calculate the terminal value (discomfort)
% Adapted from MultiInv by BryanP
% Author: Li Wang

coef_cost = [-0.477];%% Penalty factor of HVAC,WH,EV (discomfort coefficient)
% Discomfort function is calculated as the difference of the actual state
% minus desired state squared multiplied by weighting coefficients
% (coef_cost)
% For EV, if the actual SOC is greater than the desired SOC, there is no
% penalty applied.
v = (abs(s_list - params.desired_value))*coef_cost;
end

