function [p_vec, c_vec, params] = InvTransProb(old_inventory, order, t, ...
                                                future_values, params)
% INVTRANSPROB DP state transition probability, inventory problem (modular DP)
%
% as described in DP for fTransProb
%           Computes the transition probability and costs for a decision in
%           timestep (or solution step) t.
%        p_vec = a vector of probabilies such that p_vec(i) corresponds 
%                to the probability of transitioning to state s+1, given by 
%                state_iterator(i) from state s for the given decision.
%                  length(p_vec) = state_iterator('end')
%        c_vec = the cost (use negative for gain) associated with the
%                corresponding transition. In the inventory problem this is
%                the negative income from sales.
%
% Transition probabilites are based on:
%   params.p_demand    Probability of demand, as a row vector corresponding
%                      to prob of 0:length-1 demand
%
% Inventory params:
%   used here:
%      max_inv      the maximum inventory
%      p_demand     see above
%      sales_price  per unit sales price
%   not used here:
%      order_cost, unit_cost, hold_cost, term_unit_value
%
% see also:
%   DP, InvState, InvTermValue, InvCost, InvDecision, InvInit
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-10 15:47  BryanP      adapted from inventory_dp version 3
%   2  2010-05-14 08:00  BryanP      debugged probability vector
%   3  2010-05-14 10:00  BryanP      moved intialization code to InvInit()
%   4  2010-08-07 11:50  BryanP      adapted for DP v12 with:
%                                      - param returns
%                                      - future value functions

%-- initialize our output vectors
%Note: consider using sparse() here for some problem types. In this case,
%we have enough non-zeros that the full matrix is faster (at least at
%max_inv = 15
p_vec = zeros(params.max_inv+1,1);
c_vec = zeros(params.max_inv+1,1);

%-- first compute the total we have in stock
stock = old_inventory + order;

%-- determine possible sales
% first find max demand with a non-zero probability.
max_demand = length(params.p_demand)-1;

% and as a result, what are all the possible demand values
%Note: Demands is vectorized
demands = (0:max_demand)';

%match the demands to the future states
next_s = stock - demands;

% since we have no back orders, any negative states should be lumped into
% the 0 state for the probabilities, this means adding all of the probabilities greater
%than the demand we can meet, resulting in a zero or negative next state

% state0_idx = find(next_s<=0,1,'first');
% p_vec(1) = sum(params.p_demand(state0_idx:end));
p_vec(1) = sum(params.p_demand(next_s<=0));

%now we can remove the invalid (negative) future states
demands(next_s<0)=[];
next_s(next_s<0)=[];

%which allows us easily build up the contribution vector, using next_s as
%the index
c_vec(next_s + 1) = params.sales_price*demands;

%for the probability vector, we also remove the next_s=0 demand, since its 
%probability was computed above
demands(next_s==0)=[];
next_s(next_s==0)=[];

% now we can fill in the transition probability vector using the next_s
% indexes and only including those demands that we have not already counted 
% in the next_s = 0 case above.
p_vec(next_s(next_s>0)+1) = params.p_demand(1:length(demands));