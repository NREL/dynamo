function [ cost, parameters ] = InventoryCostToGoFn( current_inventory, ...
                                                  number_of_orders, ...
                                                  time, ...
                                                  parameters )
%INVENTORYCOSTTOGOFN Returns cost to go value submitting number_of_orders
% given the current_inventory

% demands = (0:this_max_sales)';
% incomes = income_fun(demands);
% next_s = s + feasible_orders(ord_idx) - demands;
% 
% %value for this (state,order) = E[income + future value - costs]
% % possible values
% pos_values = incomes + (1-disc_rate)*Values(next_s+1,t+1) - o_c(ord_idx) - h_c(ord_idx);
% % expected value, lumping any unmet demand into max sales bin
% exp_value = p_demand(:,1:this_max_sales) * pos_values(1:this_max_sales,:) + ...
%             pos_values(this_max_sales+1) ...
%                         * sum(p_demand((this_max_sales+1):length(p_demand)));

cost = 0;

end

