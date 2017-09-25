function [ number_of_units, ...
            parameters ] = InventoryStateIterator( state_index, ... 
                                                    time, ...
                                                    parameters )
%INVENTORYSTATEITERATOR Returns number of items in the inventory for a
%given state index

if state_index==0
    number_of_units = parameters.maximum_inventory+1;
    return
end

%The ith state corresponds to an inventory level of i-1, such that
%state(1)=0.
number_of_units = state_index-1;
end

