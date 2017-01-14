function [ decisions, parameters ] = InventoryDecisionSet( current_inventory, ...
                                                        time_index, ...
                                                        parameters)
%INVENTORYDECISIONITERATOR Returns post decision state
%   if decision_index = 0, then return number of decisions

number_of_units = current_inventory;
decisions = ( 0 : parameters.maximum_inventory - number_of_units );

decisions = decisions(decisions <= parameters.maximum_inventory);
parameters.decisions = parameters;
end

