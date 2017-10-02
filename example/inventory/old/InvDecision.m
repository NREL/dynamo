function [d, params] = InvDecision(s, i, t, future_values, params)
% INVDECISION DP decision iterator function for inventory problem (modular DP)
%
% as described in DP for fDecisionIterator...
%           returns a decision, represented by a variable of any type, that 
%           represents one of the possible decisions to be made at state s. 
%                i = '0 returns the number of possible decisions for the 
%                current state. Otherwize the range of i is 1:num decisions
%                for this state
%
% In this case d simply equals params.max_inv - s + i
%
% Inventory params:
%   used here:
%      max_inv   the maximum inventory
%   not used here:
%      p_demand, (sales_price, order_cost, unit_cost, hold_cost, term_unit_value)
%
% see also:
%   DP, InvState, InvTermValue, InvCost, InvTransProb, InvInit
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-10 15:47  BryanP      adapted from inventory_dp version 3
%   2  2010-05-14 20:30  BryanP      now i==0 for number of decisions
%   3  2010-08-07 11:50  BryanP      adapted for DP v12 with:
%                                      - param returns
%                                      - future value functions

if i == 0
    % add 1, because can also order zero
    d = params.max_inv - s + 1;
    return
end

%subtract one because i starts at 1
%TODO: Should be equivalent to simply i-1
d = params.max_inv - s - (i - 1);
