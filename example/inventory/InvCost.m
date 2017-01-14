function [c, params] = InvCost(old_inventory, order, t, params) 
% INVCOST DP decision cost for inventory problem (modular DP)
%
% as described in DP for fContribution 
%           the cost this period of making decision in given state, s at 
%           time/solution step t. In the inventory problem, this combines 
%           ordering & holding costs. Note that costs are negative
%
% Inventory params...
%   used here:
%      (order_cost)  fixed cost per order (default=4)
%      (unit_cost)   per_unit order cost (default=2)
%      (hold_cost)   per_unit_per period holding cost (default=1)
%   not used here:
%      max_inv, p_demand, (sales_price, term_unit_value)
%
% Defaults from Putterman 3.2.2 & 4.6.1
%
% see also:
%   DP, InvState, InvTermValue, InvTransProb, InvDecision, InvInit
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-10 15:47  BryanP      adapted from inventory_dp version 3
%   2  2010-05-11 21:17  BryanP      corrected hold cost: include inventory
%   3  2010-05-14 10:00  BryanP      moved intialization code to InvInit()
%   4  2010-08-07 11:50  BryanP      adapted for DP v12 with param returns

% ------ Compute costs ------
% Note the boolean return of order>0 mathematically evaluates to 0 for
% false and 1 for true. Hence we only are charged the order cost if we
% actually have a non-zero order

%Note: cost = order_cost + hold_cost;
c = params.order_cost*(order>0) + params.unit_cost*order ...
    + params.hold_cost*(order + old_inventory);

