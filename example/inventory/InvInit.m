function params = InvInit(max_inv, p_demand, params)
% INVINIT initialize defaults for inventory problem (modular DP)
%
% this function sets up any undefined parameters to their default values
% from Putterman 3.2.2 & 4.6.1
%
% Parameters
%      max_inv     The maximum inventory that can be held
%      p_demand    Probability of demand, as either
%                   + a row vector corresponding to prob of 0:length-1 demand
%                   + a scalar, interpreted as a poisson lambda value
%      params      An optional structure that predefines some or all of the
%                  inventory parameters to non-default values
%
% Inventory params initialized here (unless already defined in params):
%      order_cost       fixed cost per order (default=4)
%      unit_cost        per_unit order cost (default=2)
%      hold_cost        per_unit_per period holding cost (default=1)
%      sales_price      per unit sales price (default=8)
%      term_unit_value  value per unit in terminal period (default=0)
%
% see also:
%   InvCost, InvState, InvTermValue, InvTransProb, InvDecision
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-10 15:47  BryanP      adapted from inventory_dp version 3
%   2  2010-05-14 23:00  BryanP      expanded comments, added max_inv & p_demand

%Start by adding the max inventory
% Note: This has the added bonus of creating the params structure if it
% does not already exist, which saves us from explicitly checking
params.max_inv = max_inv;

%Transition Probability Defaults
if length(p_demand) == 1
    % option #1 (Mort's approach), normalize vector to 1
    %  	p_demand = poisspdf(0:max_inv,p_demand);
    %  	p_demand= p_demand/sum(p_demand);

    % option #2 (Bryan's approach), lump any additional demand into the
    % highest demand... this seems consistant with the idea that any
    % additional demand can only produce the same ammount of sales as the
    % maximum inventory
 	params.p_demand = poisspdf(0:(max_inv-1),p_demand)';
 	params.p_demand(max_inv+1) = 1 - sum(params.p_demand);
else
    params.p_demand = p_demand;
end


if not(isfield(params, 'sales_price'))
    params.sales_price = 8;
end

%Terminal Value Defaults
if not(isfield(params, 'term_unit_val'))
    params.term_unit_val = 0;
end

%Cost Defaults
if not(isfield(params, 'order_cost'))
    params.order_cost = -4;
end

if not(isfield(params, 'unit_cost'))
    params.unit_cost = -2;
end

if not(isfield(params, 'hold_cost'))
    params.hold_cost = -1;
end

