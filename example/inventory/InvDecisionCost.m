function c_list = InvDecisionCost(params, t, pre_state_list, order_list)  %#ok<INUSL> extra params listed by name for clarity
% INVDECISIONCOST DP decision cost for inventory problem
%
% Usage: c_list = InvDecisionCost(params, t, pre_state_list, order_list) 
%
%           the cost this period of making decision in given state, s at 
%           time/solution step t. In the inventory problem, this includes a
%           per order cost (if any) plus the sum of per item costs.
%           Note that costs are negative. 
%
% Multiple state/orders are supported by stacking state and order entries
% to form column vectors with one state/order per row
%
% see also:
%
% adapted from MultiInvDecisionCost v10 by Bryan Palmintier 2017

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-09-24 21:09  BryanP      adapted from MultiInvDecisionCost v10 

% ------ Compute costs ------
% Note the boolean return of sum(order)>0 mathematically evaluates to 0 for
% false and 1 for true. Hence we only are charged the order cost if we
% actually have a non-zero order

%Note: cost = order_cost + hold_cost;
c_list = params.order_cost*(sum(order_list,2)>0) ...
    + order_list * params.unit_cost';

end

