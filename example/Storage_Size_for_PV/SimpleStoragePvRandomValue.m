function c_list = SimpleStoragePvRandomValue(params, t, post_state_list, uncertain_demand_list) %#ok<INUSD,INUSL>
% SIMPLESTORAGEPVRANDOMVALUE DP uncertainty cost for simple sotrage + PV problem
%
% Usage: c_list = SimpleStoragePvRandomCost(params, t, post_state_list, uncertain_demand_list)
%
% There are no costs to the utility for the random process in this problem.
%
% Note: each state or order entry must be a row, and both should have
% params.n_products columns.
%
% see also:
%
% originally by Bryan Palmintier 2016

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-07-18 07:06  BryanP      adapted from MultiInvRandomCost v2

% ------ Compute costs ------
c_list = zeros(size(max(post_state_list,uncertain_demand_list), 1), 1);

