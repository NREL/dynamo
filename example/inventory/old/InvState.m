function [s, params] = InvState(i, t, params)
% INVSTATE DP state iterator function for inventory problem (modular DP)
%
% as described in DP for fStateIterator...
%           returns a state variable (of any type), for the positive 
%           decimal state index "i"
%                i = 0 returns the number of possible states
%                otherwise the range of i is 1:number of states
%
% In this case s simply equals i-1; and InvState(0,params) = max_inv+1 (for zero)
%
% Inventory params:
%   used here:
%      max_inv   the maximum inventory, 
%   not used here:
%      p_demand, (sales_price, order_cost, unit_cost, hold_cost, term_unit_value)
%
% see also:
%   DP, InvDecision, InvTermValue, InvCost, InvTransProb, InvInit
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-10 15:47  BryanP      adapted from inventory_dp version 3
%   2  2010-05-14 20:30  BryanP      now i==0 for number of states
%   3  2010-05-22 09:50  BryanP      added support for time-aware states
%   4  2010-08-07 11:50  BryanP      adapted for DP v12 with param returns


% When called with the state index of 'end' we need to return the total
% number of states. Which in this case is max_inv+1, since we also have a
% zero state.
if i==0
    s = params.max_inv+1;
    return
end

%The ith state corresponds to an inventory level of i-1, such that
%state(1)=0.
s = i-1;
