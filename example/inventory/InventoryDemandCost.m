function [cost, parameters] = InventoryDemandCost( state, decision, period, values, parameters )
%INVENTORYDEMANDCOST Summary of this function goes here
%   Detailed explanation goes here

%-- initialize our output vectors
%Note: consider using sparse() here for some problem types. In this case,
%we have enough non-zeros that the full matrix is faster (at least at
%max_inv = 15
p_vec = zeros(parameters.maximum_inventory+1,1);
c_vec = zeros(parameters.maximum_inventory+1,1);


cost = 0;

end

