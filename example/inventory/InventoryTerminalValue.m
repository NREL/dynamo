function [ terminal_value, parameters ] = InventoryTerminalValue( state_index, ...
                                                            parameters)
number_of_units = state_index-1;
terminal_value = number_of_units * parameters.terminal_unit_value;
end

