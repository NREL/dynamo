function [ cost, parameters ] = InventoryDecisionCost( current_inventory, number_of_orders, current_period, parameters )

cost = parameters.order_cost * (number_of_orders);

end

