function c_list = MultiInvDecisionCost(params, t, pre_state_list, order_list)  %#ok<INUSL> extra params listed by name for clarity
% MULTIINVDECISIONCOST DP decision cost for multi-product inventory problem
%
% Usage: c_list = MultiInvDecisionCost(params, t, pre_state_list, order_list) 
%
%           the cost this period of making decision in given state, s at 
%           time/solution step t. In the inventory problem, this includes a
%           per order cost (if any) plus the sum of per item costs.
%           Note that costs are negative. 
%
% Note: each state or order entry must be a row, and both should have
% params.n_products columns.
%
% Multiple state/orders are supported by stacking state and order entries
% (each a row)
%
% see also:
%
% originally by Bryan Palmintier 2010
% reworked for shared problem structure by Bryan Palmintier 2016

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  10  2016-10-27 14:41  BryanP      Removed params return, corrected usage 
%   9  2016-10-27 12:31  BryanP      Reordered parameters to standardize on t as second input 
%   8  2016-10-27 12:05  BryanP      Moved hold cost to MultiInvOps 
%   7  2016-10-21 15:40  BryanP      Reordered input parameters to put params first (for possible future objectifying) 
%   6  2016-07-07 00:30  BryanP      Renamed MultiInvCost to MultiInvDecisionCost for shared problem structure 
%   5  2016-05-01 04:20  BryanP      Allow one state multi-decision, and one decision, mulit-state through bsxfun 
%   4  2016-04-29 21:40  BryanP      Convert to row vector for each state and support vectors of states 
%   3  2010-08-08 00:30  BryanP      adapted for DP v12 with param returns
%   2  2010-05-15 20:45  BryanP      updated to use params.states
%   1  2010-05-15 09:47  BryanP      adapted from InvCost v3

% ------ Compute costs ------
% Note the boolean return of sum(order)>0 mathematically evaluates to 0 for
% false and 1 for true. Hence we only are charged the order cost if we
% actually have a non-zero order

%Note: cost = order_cost + hold_cost;
c_list = params.order_cost*(sum(order_list,2)>0) ...
    + order_list * params.unit_cost';

