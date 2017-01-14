function [Orders, Values] = inventory1(max_inv, n_periods, p_demand)
% INVENTORY1 Simple inventory (ESD.862 HW4 #1a)
%
% Usage: [Orders, Values] = inventory_dp(max_inv, n_periods, p_demand)
%
%   Required
%     max_inv     Maximum Inventory (scalar)
%     n_periods   Number of time periods to consider (scalar)
%   Optional
%     p_demand    Probability of demand, as either
%                   + a row vector corresponding to prob of 0:length-1 demand
%                   + a scalar, interpreted as a poisson lambda value
%     default functions & p_demand from Putterman 3.2.2 & 4.6.1
%
% Solves the stochastic inventory problem using traditional finite horizon
% dynamic programming via backwards induction.
%
% See also inventory_dp
%
% Originally by Bryan Palmintier, 2010


% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-11-01 11:50  BryanP      adapted from inventory_dp ver 3

if nargin <3
    p_demand = [0.25, 0.5, 0.25];
end

[Orders, Values] = inventory_dp(max_inv, n_periods, p_demand);
