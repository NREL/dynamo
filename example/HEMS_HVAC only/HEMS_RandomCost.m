function [ c_list ] = HEMS_RandomCost( params, t, post_state_list, uncertain_demand_list )
%function [ c_list ] = HEMS_RandomCost( params, t, post_state_list, uncertain_demand_list )
% Calculate the random cost for EV only. HVAC and WH do not have random
% cost.
% Adapted from MultiInv by BryanP
% Author: Li Wang

%only HVAC
c_list = zeros(size(uncertain_demand_list(:,1)));

% 
% coef_ele = params.electricity_price(t);  % Price of electricity, $/KWh
% EV_l_bound = params.appliance_range{3}(1,1); % Lower bound for EV SOC
% 
% %If the EV is at home and the post-decision state does not hit the lower
% %bound of EV SOC, there is no random cost. Otherwise, there is electriciy
% %cost for non-controllable load to be supplied by grid directly.
% if (post_state_list(3) > EV_l_bound && params.arrive_departure_time(1,1) < t && t <= params.arrive_departure_time(1,2))
%     c_list = zeros(size(uncertain_demand_list(:,1))); % No random cost when the EV is at home and does not hit the lower bound.
% else
%     c_list = uncertain_demand_list(:,3)*params.unit_cost(3)*coef_ele; % Cost for supplying non-controllable load when EV is not at home or hits the lower bound.
% end

