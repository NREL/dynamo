function [v, params] = InvTermValue(s, t, params)
% INVTERMVALUE DP terminal value for inventory problem (modular DP)
%
% as described in DP for fFinalValue
%           returns the final value for state s
%
% Inventory params...
%   used here:
%      term_unit_value   value per unit in terminal period
%   not used here:
%      max_inv, p_demand, sales_price, order_cost, unit_cost, hold_cost
%
% Defaults from Putterman 3.2.2 & 4.6.1
%
% see also:
%   DP, InvState, InvCost, InvTransProb, InvDecision, InvInit
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-10 15:47  BryanP      adapted from inventory_dp version 3
%   2  2010-05-14 10:00  BryanP      moved intialization code to InvInit()
%   3  2010-08-07 11:50  BryanP      adapted for DP v12 with param returns

% ------ Compute costs ------
v = s * params.term_unit_val;
