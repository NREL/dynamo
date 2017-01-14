function [p_vec, c_vec, params, next_s_list] = MultiInvTransProb(s, order, t, ...
                                    future_values, params) %#ok<INUSL>
% MULTIINVTRANSPROB DP state transition probability for multi-product multi-product inventory problem 
%
% as described in DP for fTransProb
%           Computes the transition probability and costs for a decision in
%           timestep (or solution step) t.
%        p_vec = a vector of probabilies such that p_vec(i) corresponds 
%                to the probability of transitioning to state s+1, given by 
%                state_iterator(i) from state s for the given decision.
%                  length(p_vec) = state_iterator('end')
%        c_vec = the cost (use negative for gain) associated with the
%                corresponding transition. In the multi-product multi-product inventory problem this is
%                the negative income from sales.
%        Note: t and future_values are not used
%
% see also:
%   DP, MultiInvState, MultiInvTermValue, MultiInvCost, MultiInvDecision, MultiInvInit
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-15 09:47  BryanP      adapted from InvTransProb v3
%   2  2010-06-07 00:30  BryanP      Allow param update with each sub-function call 
%   3  2011-03-29 10:10  BryanP      Renamed MultiInvValidCombin to CombinWithLimits
%   4  2016-04-14 11:23  BryanP      Enable returning corresponding next states
%   5  2016-04-15 02:03  BryanP      Enable optional export of next states  
%   6  2016-04-29 22:20  BryanP      Convert to row vector for each state 

%-- initialize our output vectors
%Note: consider using sparse() here for some problem types. In this case,
%we have enough non-zeros that the full matrix is faster.
p_vec = zeros(params.n_states,1);
c_vec = zeros(params.n_states,1);

%-- rename the state
old_inventory = s;

%-- first compute the total we have in stock
stock = old_inventory + order;

% -- Now compute the possible future states

% Rather than listing demands and computing states, here we opt to work 
% backwards and list the future states and back solve for the demands that
% got us there. 
next_s_list = PossibleFuture(stock, [params.demand(:).max]);

% Compute corresponding demands
demands = ones(size(next_s_list,1),1)*stock - next_s_list;

%-- Compute the income from these sales
c_vec = demands * params.sales_price';

%-- Find probabilities for each demand level
%
% Future note: for n-dimensional joint probabilities, use sub2ind() to
% handle all cases with demand <= inventory and then identify the subset of
% demands requiring special treatment using min().
p_vec = ones(size(c_vec));  %start as ones so can multiply to find joint probabilities 

%- first compute the marginal probabilities required to fully sell out
%  of each item (separately) by summing any possible demands that exceed the
% current inventory
p_sellout = zeros(1, params.num_products);

for p = 1:params.num_products;
    %Sum all probabily of selling >= the current inventory--indexed by
    %stock(p) + 1, because prob(1) is for zero demand.
    p_sellout(p) = sum(params.demand(p).prob( (stock(p)+1) : (params.demand(p).max+1)));
end

%- next match marginal (per product) probabilities for non-sellout demands
for p = 1:params.num_products;
    % Vectorize by masking for all non-sellout demands. Included for
    % clarity (strictly speaking this is may not be needed since non-zeros
    % evalutate to true)
    sellout_mask = next_s_list(:,p) == 0;
    
    %Use previously calculated sellout probabilitY (scalar) where needed
    p_vec(sellout_mask) = p_vec(sellout_mask) * p_sellout(p);
    % Otherwise, use the corresponding column of the demand as indices into
    % the probability lookup table and update the corresponding entries.
    % Note must again increment the actual demands by one since first
    % element of prob corresponds to zero demand
    if any(not(sellout_mask))
        p_vec(not(sellout_mask)) = p_vec(not(sellout_mask)) ...
                    .* params.demand(p).prob( demands(not(sellout_mask), p) + 1);
    end
end


% %OLD
% %- now compute the joint probability (assuming independance) for each
% %demand level
% for d = 1:size(demands,1);
%     this_demand = demands(d,:);
%     prob = 1;
%     for product = 1:params.num_products;
%         prob = prob  * p_sales{product}(this_demand(product)+1);
%     end
%     
%     p_vec(d) = prob;
% end
assert(abs(sum(p_vec) - 1)<1e-6, 'Joint probability does not sum to 1')

%Convert next state list to a cell vector with one row vector per cell.
next_s_list = mat2cell(next_s_list, ones(size(next_s_list,1),1));

    
end %main function

function next_s = PossibleFuture(stock, max_demand)
% Recursive helper function to list possible future states
%
% Usage: next_s = PossibleFuture(stock, max_demands)
%
% Inspired by CombinWithLimits() v. 1

    %Basecase (ends the recusion)
    if length(stock) == 1;
        next_s = (max(0,stock-max_demand):stock)';
        return;
    end
    
    %initialize the output array
    next_s = [];
    
    %Handle the first item posibilities here, and other posibilities via
    %recursion
    %
    % For our first item, the possible future states will run from zero
    % to our current stock, unless the maximum demand is not high enough,
    % in which case, we would only run from the lowest inventory level we
    % could achieve
    for num_this_item = max(0,stock(1)-max_demand(1)):stock(1)
        
        %recursively get sub-combinations
        other_items = PossibleFuture(stock(2:end),max_demand(2:end));
        new_rows = size(other_items,1);
        
        %now add these new combinations to our list of item inventories
        %Note: that we prepend the num_this_item to the new other_items,
        %such that the top level state list will have an entry for every
        %product
        next_s = [ next_s
                   [ones(new_rows,1)*num_this_item  other_items]
                 ]; %#ok<AGROW>
    end
end

