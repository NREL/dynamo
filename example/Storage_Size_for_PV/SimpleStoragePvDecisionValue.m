function value_list = SimpleStoragePvDecisionValue(params, t, pre_state_list, decision_list)
% SIMPLESTORAGEPVDECISIONCOST DP decision cost for simple storage plus PV problem
%
% Usage: value_list = SimpleStoragePvDecisionCost(params, t, pre_state_list, decision_list) 
%
%           the cost this period of making decision in given state, s at 
%           time/solution step t. 
%
% Nominally this is the inverter plus storage costs for the new capacity,
% But we also want to enforce unit sizes for each, so we first check to see
% if we have a unit size match and adjust accordingly
% 
% Notes:
%  -- it is possible to scrap a system (decrease size) in which case the
%     value is heavily discounted
%  -- costs have negative value
%
% Note: each state and order entry must be a row
%
% Multiple state/orders are supported by stacking state and order entries
% (each a row)
%
% see also:
%
% originally by Bryan Palmintier 2017

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-07-16 22:37  BryanP      adapted from MultiInvCost v10

% Ensure state and decision list are same length
n_decision = size(decision_list, 1);
n_state = size(pre_state_list, 1);

% Duplicates states accordingly
if n_state < n_decision
    pre_state_list = repmat(pre_state_list, n_decision, 1);
elseif n_decision < n_state
    decision_list = repmat(decision_list, n_state, 1);
end

% ------ Compute costs ------
% Determine system size to change based on need for replacement or not
new_build_list = decision_list;
% ID inverter size increase and adjust new_build if require full
% system replacement (ISGT assumption)
inv_upsize_map = (decision_list(:,1) > 0);
new_build_list(inv_upsize_map,:) = ...
    bsxfun(@times, new_build_list(inv_upsize_map, :) + pre_state_list(inv_upsize_map, 2:3), params.replace_on_inv_upsize);
% Finally remove any scrapping from the new build list (b/c costed
% different)
new_build_list = max(0, new_build_list);

%ID scrap quantities
scrap_list = zeros(size(decision_list));
scrap_map = (decision_list < 0);
scrap_list(scrap_map) = abs(decision_list(scrap_map));
%Adjust numbers when working with full replacements
inv_downsize_map = (decision_list(:,1) < 0);
scrap_list(inv_downsize_map,:)  = bsxfun(@times, abs(pre_state_list(inv_downsize_map,2:3)), params.scrap_value_on_inv_upsize);

%--Finally, compute the corresponding values
% Start with new + scrap% * scrap
feature_cost = [params.inv_cost(t) + params.bos_cost(t), params.bat_cost(t)];
value_list = bsxfun(@times, new_build_list, feature_cost) ...
                - params.scrap_ratio * bsxfun(@times, scrap_list, feature_cost);
% Compute transformer values
%  Note: no scrap value for transformer (for simplicity)
transformer_values =  params.xfmr_cost_base(t) * max(0, new_build_list(:,1)) .^ params.xfmr_cost_expon(t);

value_list = sum(value_list, 2) + transformer_values;