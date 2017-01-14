function [Orders, Values] = inventory2(max_inv, n_periods, p_demand)
% INVENTORY2 Simple inventory with Poisson demand(ESD.862 HW4 #1b)
%
% Usage: [Orders, Values] = inventory_dp(max_inv, n_periods, p_demand)
%
%   Required
%     max_inv     Maximum Inventory (scalar)
%     n_periods   Number of time periods to consider (scalar)
%     p_demand    Probability of demand, as either
%                   + a row vector corresponding to prob of 0:length-1 demand
%                   + a scalar, interpreted as a poisson lambda value
%     default functions & p_demand from Putterman 3.2.2 & 4.6.1
%
% Solves the stochastic inventory problem using traditional finite horizon
% dynamic programming via backwards induction.
%
% See also inventory_dp, inventory1
%
% Originally by Bryan Palmintier, 2010


% Implementation Note: Yes, this is identical to inventory1, but without the
% default p_demand...
%
% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-11-01 11:50  BryanP      adapted from inventory_dp ver 3

[Orders, Values] = inventory_dp(max_inv, n_periods, p_demand);
