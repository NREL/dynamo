function c_list = InvRandomCost(params, t, post_state_list, uncertain_demand_list) %#ok<INUSL>
% INVRANDOMCOST DP uncertainty cost for multi-product inventory problem
%
% Usage: c_list = InvRandomCost(params, t, post_state_list, uncertain_demand_list)
%
% Here uncertainty is in demand and hence money earned. This is captured as
% a negative cost (as captured in the params)
%
% Note: each state or order entry must be a row, and both should have
% params.n_products columns.
%
% Multiple state/orders are supported by stacking state and order entries
% (each a row)
%
% see also:
%
% originally by Bryan Palmintier 2016

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   2  2016-12-01 10:26  BryanP      BUGFIX: use bsxfun to handle different sized lists 
%   1  2016-10-27 12:31  BryanP      adapted from InvDecisionApply v9

% ------ Compute costs ------
%Note: have to cap the contribution to only consider fulfilled demand by
%capping at current (inventory) state if needed
c_list = bsxfun(@min, post_state_list, uncertain_demand_list) * params.sales_price';

