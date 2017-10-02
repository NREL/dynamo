clc
clear

params.order_cost = 1;
params.terminal_unit_value = 1;
params.maximum_inventory = 4;

problem.n_periods = 3;
problem.state_set = {(1:4)',(1:4)',(1:4)'};
problem.fTerminalValue = @InventoryTerminalValue;
problem.fDecisionSet = @InventoryDecisionSet;
problem.fDecisionApply = @InventoryDecisionCost;
problem.fUncertaintyApply = @InventoryDemandCost;
problem.params = params;

results = dpBI(problem);           
