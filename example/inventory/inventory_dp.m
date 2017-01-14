function [Orders, Values, p_demand] = inventory_dp(max_inv_or_struct, n_periods, p_demand, disc_rate, ...
								order_cost_fun, income_fun, hold_cost_fun, term_reward_fun)
% INVENTORY_DP Solve Dynamic Programing Stochastic Inventory Problem (ESD.862)
%
% Call Options: 
%   INV_STRUCT MODE: [Orders, Values, p_demand] = inventory_dp(inv_prob_struct, n_periods)
%   ESD.862 Mode:    [Orders, Values, p_demand] = inventory_dp(max_inv, n_periods, p_demand, disc_rate ...
%								@order_cost_fun, @income_fun, @hold_cost_fun, @term_reward_fun)
%
%   INV_STRUCT MODE, 
%   all of the problem parameters are provided in the inv_prob_struct and
%   any missing non-required fields (as defined in ESD.862 mode) revert to
%   the ESD.862 defauts. This provides compatibility with inventory_adp
%   IMPORTANT in inv_struct mode, field names are deterimined by
%   inventory_adp:
%       inv_prob.lambda      Specifies poisson arrival rate (prob vector
%                            not supported)
%       inv_prob.dr          Discount rate [optional]
%       inv_prob.max_inv     Maximum Inventory (scalar)
%       inv_prob.n_periods   Number of time periods to consider (scalar),
%                            optionally specified/overridden with n_periods
%                            parameter to inventory_dp.
%   And, notably, the cost function parameters are defined rather than the
%   cost functions themselves. These are:
%       inv_prob.order_cost_fix  Fixed cost per order [4/order]
%       inv_prob.order_cost_ea   Per unit ordering cost [2/ea]
%       inv_prob.sales_price     Income on sales [8/ea]
%       inv_prob.hold_cost       Cost to keep in inventory [1/ea]
%       inv_prob.term_reward     Value of any inventory in final period [0/ea]
%
%
%   ESD.862 MODE, 
%     max_inv     Maximum Inventory (scalar)
%     n_periods   Number of time periods to consider (scalar)
%     p_demand    Probability of demand, as either
%                   + a row vector corresponding to prob of 0:length-1 demand
%                   + a scalar, interpreted as a poisson lambda value
%   Optional
%     disc_rate   The discount rate per period. (default = 0)
%     
%     functions for cost(units), income(units sold), hold(units held), terminal reward(units)
%     default functions from Putterman 3.2.2 & 4.6.1
%
% Solves the stochastic inventory problem using traditional finite horizon
% dynamic programming via backwards induction.
%
% Note: The function combines the inventory1, inventory2, and inventory3
% functions defined in ESD.862 HW4
%
% See also inventory1, inventory2, inventory3, inventory_adp
%
% Originally by Bryan Palmintier, 2009


% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2009-Fall         BryanP      Practice problems for generals
%   2  2009-Fall         BryanP      added support for poisson demand
%   3  2010-04-27 22:17  BryanP      added discounting & extra comments
%   4  2010-11-03 22:40  BryanP      Added support for ADP-like parameter structures

if isstruct(max_inv_or_struct)
    inv_prob = max_inv_or_struct;
    
    %--- Extract parameters
    % first the required ones
    max_inv = inv_prob.max_inv;
    p_demand = inv_prob.lambda;
    % use n_periods parameter if given
    if nargin == 1
        n_periods = inv_prob.n_periods;
    end
    
    % then the optional ones
    if isfield(inv_prob, 'dr')
        disc_rate = inv_prob.dr;
    end
    
    %and finaly cost parameters to conver to cost functions
    if isfield(inv_prob, 'sales_price')
    	income_fun = @(sold) inv_prob.sales_price * sold;
    end
    if isfield(inv_prob, 'hold_cost')
    	hold_cost_fun = @(held) inv_prob.hold_cost * held;
    end
    if isfield(inv_prob, 'term_reward')
    	term_reward_fun = @ (units) inv_prob.term_reward *units;
    end
    
    % with special care for the two parameter order cost function
    % Such that any combination of non-defaults can be specified
    if not(isfield(inv_prob, 'order_cost_fix'))
    	inv_prob.order_cost_fix = 4;
    end
    if not(isfield(inv_prob, 'order_cost_ea'))
    	inv_prob.order_cost_ea = 2;
    end
	order_cost_fun = @(units) inv_prob.order_cost_fix * (units>0) ...
                                + inv_prob.order_cost_ea * units;
    
else
    max_inv = max_inv_or_struct;
end

% --- Setup Problem ---
% Setup poisson demand probability vector if needed
if length(p_demand) == 1
    % option #1 (Mort's approach), normalize vector to 1
    %  	p_demand = poisspdf(0:max_inv,p_demand);
    %  	p_demand= p_demand/sum(p_demand);

    % option #2 (Bryan's approach), lump any additional demand into the
    % highest demand... this seems consistant with the idea that any
    % additional demand can only produce the same ammount of sales as the
    % maximum inventory
 	p_demand = poisspdf(0:(max_inv-1),p_demand);
 	p_demand(max_inv+1) = 1 - sum(p_demand);
end

% -- Default to no discounting
if nargin < 4 || not(exist('disc_rate', 'var')) || isempty(disc_rate)
	disc_rate = 0;
end

% -- Default functions
%ordering cost
if nargin < 5 || not(exist('order_cost_fun', 'var')) || isempty(order_cost_fun)
	order_cost_fun = @(units) 4*(units>0) + 2*units;
end
%income from unit sales
if nargin < 6 || not(exist('income_fun', 'var')) || isempty(income_fun) 
	income_fun = @(sold) 8*sold;
end

%cost of inventory held between periods
if nargin < 7 || not(exist('hold_cost_fun', 'var')) || isempty(hold_cost_fun) 
	hold_cost_fun = @(held) 1*held;
end

%terminal reward function
if nargin < 8 || not(exist('term_reward_fun', 'var')) || isempty(term_reward_fun)
	term_reward_fun = @ (units) 0*units;
end

% -- Parameter intialization
States = (0:max_inv)'; %vector of possible states, index = inventory+1
n_states = length(States);
max_demand = (length(p_demand)-1);

%value function = f(s,t). Initialize to -Inf so any real value is greater
Values = -Inf*ones(n_states, n_periods+1);
Orders = NaN * ones(n_states, n_periods); %action (orders) = f(s,t), no action in final state
	

% --- Work backwards ---
% -- Start with terminal period
Values(:,n_periods+1) = term_reward_fun(States);

% Loop over time periods in reverse (backward induction)
for t = n_periods:-1:1
    % Loop over states
	for s_idx = 1:n_states;
		s=States(s_idx);
		
		%Note: feasible orders is vectorized (as long as possible)
        %inventory capacity limits the feasible orders
		feasible_orders = 0:(max_inv-s);
        %compute resulting costs (order & hold) for all feasible orders
		o_c = order_cost_fun(feasible_orders);
		h_c = hold_cost_fun(s+feasible_orders);
		%compute maximum possible sales for each order amount
		max_sales_vec = min(s+feasible_orders, max_demand); %vector
		
		%Now Loop over each possible order
		for ord_idx = 1:length(feasible_orders)
            %Note: Demand is vectorized
			this_max_sales = max_sales_vec(ord_idx);
			demands = (0:this_max_sales)';
			incomes = income_fun(demands);
			next_s = s + feasible_orders(ord_idx) - demands;
			
			%value for this (state,order) = E[income + future value - costs]
			% possible values
			pos_values = incomes + (1-disc_rate)*Values(next_s+1,t+1) - o_c(ord_idx) - h_c(ord_idx);
			% expected value, lumping any unmet demand into max sales bin
			exp_value = p_demand(:,1:this_max_sales) * pos_values(1:this_max_sales,:) + ...
						pos_values(this_max_sales+1) ...
									* sum(p_demand((this_max_sales+1):length(p_demand)));
			if  exp_value > Values(s+1, t)
				Values(s_idx, t) = exp_value;
				Orders(s_idx, t) = feasible_orders(ord_idx);
			end
		end
		
	end 	%End (loop over all states s for time t)

end % End (loop over all time periods)
