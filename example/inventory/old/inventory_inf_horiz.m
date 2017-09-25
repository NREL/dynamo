function [Orders, Values, n_iter, v_trace] = inventory_inf_horiz(max_inv, p_demand, ...
        disc_rate, type, tol, guess, order_cost_fun, income_fun, hold_cost_fun)
% INVENTORY_INF_HORIZ Infinite Horizon Stochastic Inventory DP (ESD.862)
%
% Usage: [Orders, Values, n_iter, v_trace] = inventory_inf_horiz(max_inv, p_demand,  ...
%        disc_rate, type, tol, guess, order_cost_fun, income_fun, hold_cost_fun)
%   Required
%     max_inv     Maximum Inventory (scalar)
%     p_demand    Probability of demand, as either
%                   + a row vector corresponding to prob of 0:length-1 demand
%                   + a scalar, interpreted as a poisson lambda value
%                   
%   Optional
%     disc_rate   The discount rate per period. (default = 0)
%     type        Type of iteration. Options include:
%                    'value'  converge the value vector(default)
%                    'policy' converge the order policy vector via markov
%                             policy iteration
%                    'val till pol' DOES NOT WORK uses value iteration, 
%                             but checks for policy convergence for 2
%                             sequential periods. trouble is DOES NOT WORK
%     tol         the (relative) convergence tolerence (default = 0.001)
%     guess       vector of starting Value for iteration. (default = all zeros)
%     order_cost_fun, income_fun, hold_cost_fun...
%                 functions for cost(units), income(units sold), 
%                 hold(units held).Defaults from Putterman 3.2.2 & 4.6.1
%
% Solves the stochastic inventory problem for infinite horizon
% using dynamic programming via either value or policy iteration.
%
% Originally by:
%   Bryan Palmintier
%   for ESD.862, Spring 2010


% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-04-28 15:40  BryanP      adapted from inventory_dp ver 3
%   2  2010-04-29 12:00  BryanP      restructured for other iteration types
%   3  2010-04-29 23:00  BryanP      Policy Iteration working... finally!
%   3  2010-04-29 23:50  BryanP      Added ability to trace value evolution

%% --- Setup Problem ---
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

% --- Handle defaults ---
% -- Default to no discounting
if nargin < 3
	disc_rate = 0;
end
% -- Setup for iteration
% iteration type
if nargin < 4
    type = 'value';
end

% tolerance
if nargin < 5
    tol = 0.001;
end

% initial guess
if nargin < 6
    guess = zeros(max_inv+1,1);
end

% -- Default functions
%ordering cost
if nargin < 7
	order_cost_fun = @(units) 4*(units>0) + 2*units;
end
%income from unit sales
if nargin < 8
	income_fun = @(sold) 8*sold;
end

%cost of inventory held between periods
if nargin < 9
	hold_cost_fun = @(held) 1*held;
end

% ---- Intialize Parameters ----
%-- program control flags
%create flag for which type of iteration to use
use_value_iteration = strcmp(type, 'value') || strcmp(type, 'val till pol');
check_policy = strcmp(type, 'val till pol') || strcmp(type, 'policy');
store_trace = nargout > 3;

%-- state space
States = (0:max_inv)'; %vector of possible states, index = inventory+1
n_states = length(States);
max_demand = (length(p_demand)-1);

%-- value & policy vectors
%save old vector for comparision
if use_value_iteration
    Orders = NaN * ones(n_states, 1); %action (orders) = f(s)
    Values = guess;
    if check_policy
        old = [Orders, Orders];
    else
        old = Values;
    end
else %policy iteration
    Orders = guess;
    Values = NaN * ones(1, n_states);
    old = Orders;
end

%-- convergence criteria
% scale tolerance for discount factor
tol = tol*disc_rate/(2*(1-disc_rate));
%initialize error
err = Inf;
n_iter = 0;
    
% >>>>>>>> The main loop <<<<<<<<<	
% Iterate until converged
while (err > tol)
    
    %manage iterations
    n_iter = n_iter+1;
    
    %Use the desired iteration method
    if use_value_iteration
        %% ------ Value Iteration ------
        % Implementation Note: this does scream for a new function, but the
        % parameter list would be annoyingly long, so here it is in-line
        
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
            for ord_idx = 1:length(feasible_orders);
                %Note: Demand is vectorized
                this_max_sales = max_sales_vec(ord_idx);
                demands = (0:this_max_sales)';
                incomes = income_fun(demands);
                next_s = s + feasible_orders(ord_idx) - demands;

                % >>>>>> Core change for infinte horizon is here <<<<<<
                % The difference is that the next future state is not the next
                % column in the value matrix, but rather refers to the same
                % state vector we are already using, b/c of steady-state
                % conditions

                % Implementation note:
                %   - This algorithm uses Gauss-Siedel updating, such that the
                %   updated value functions are immediately available, rather
                %   than waiting ubtil the next iteration for faster
                %   convergence

                %value for this (state,order) = E[income + future value - costs]
                % possible values
                pos_values = incomes + (1-disc_rate)*Values(next_s+1) - o_c(ord_idx) - h_c(ord_idx);
                % expected value, lumping any unmet demand into max sales bin
                exp_value = p_demand(:,1:this_max_sales) * pos_values(1:this_max_sales,:) + ...
                            pos_values(this_max_sales+1) ...
                                        * sum(p_demand((this_max_sales+1):length(p_demand)));
                %if this order has a better expected value, use it to update
                %the value for this state.
                if  exp_value > Values(s+1)
                    Values(s_idx) = exp_value;
                    Orders(s_idx) = feasible_orders(ord_idx);
                end
            end %Loop over orders
        end %End (loop over all states s for this iteration)

        if check_policy
            err = norm(Orders-old(:,1)) + norm(Orders-old(:,2));
            old = [Orders, old(:,1)];
            if isnan(err)
                err = Inf;
            end
        else
            err = norm(Values-old);
            old = Values;
        end

    else %policy iteration
        %% ------ Policy Iteration ------
        % Implementation Note: this does scream for a new function, but the
        % parameter list would be annoyingly long, so here it is in-line
        
        %algorithmic steps numbered to match Powell (2007) Figure 3.5 (p. 62)
        %1a> compute the one step transition matrix
        [P, ev_incomes] = trans_matrix(Orders, max_inv, 0:max_inv, p_demand, income_fun);
        
        %1b> compute the "reward" function for the state at current time
        net_income = ev_incomes - order_cost_fun(Orders) ... 
                        - hold_cost_fun(Orders + (0:max_inv)');
        
        %2> use this to estimate the new value fuction using a vector form
        %of Belman's equation with discounting. The equation is:
        %      V_n = income + r * E[V_n+1|X]
        %but then we rearrange as
        %      income = V_n - r * E[V_n+1|X]
        %             = V_n - r * P[V_n+1|X]V_n+1
        %recalling V_n = V_n+1 for infinite horizon we have
        %      income = (1-r*P[V|X])*Values
        %solving for Values we have:
        %      Values = income/(1-r*P[V|X])
        %then we convert to vector notation to solve all states
        %simultaneously...
        Values = (eye(n_states) - (1-disc_rate)*P)\net_income;
        
        %3> find the new optimal policy (Orders) given these Values
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
            best_exp_value = -Inf;
            for ord_idx = 1:length(feasible_orders);
                %Note: Demand is vectorized
                this_max_sales = max_sales_vec(ord_idx);
                demands = (0:this_max_sales)';
                incomes = income_fun(demands);
                next_s = s + feasible_orders(ord_idx) - demands;

                %value for this (state,order) = E[income + future value - costs]
                % possible values
                pos_values = incomes + (1-disc_rate)*Values(next_s+1) - o_c(ord_idx) - h_c(ord_idx);
                % expected value, lumping any unmet demand into max sales bin
                exp_value = p_demand(:,1:this_max_sales) * pos_values(1:this_max_sales,:) + ...
                            pos_values(this_max_sales+1) ...
                                        * sum(p_demand((this_max_sales+1):length(p_demand)));
                % >>>>>> Core change for policy iteration is here <<<<<<
                % Compared to value iteration, the difference here is that
                % we are ONLY saving the best policy (Order) decision 
                if  exp_value > best_exp_value
                    Orders(s_idx) = feasible_orders(ord_idx);
                    best_exp_value = exp_value;
                end
            end %Loop over orders
        end %End (loop over all states s for this iteration)
        
        
        err = norm(Orders-old);
        old = Orders;
        
    end %iteration type if
    
    if store_trace
        v_trace(:,n_iter) = Values; %#ok<AGROW>
    end

end % End (loop over all time periods)
end %Inventory_inf_horiz function

%% ----- Helper Functions -----
function [P, ev_incomes] = trans_matrix(Orders, max_inv, demand, demandprob, income_fun)
% Helper function to build the single stage policy transition matrix, P
%
%  Note: demand and demandprob must correspond
%
% i-th row is the probability of transitioning from state i (inventory
% level of i-1) to all other states (one per column) given the current
% policy for this state (Orders(i)) and the demand probabilities

    %preallocate our outputs
    P = zeros(max_inv + 1, max_inv + 1);
    ev_incomes = zeros(max_inv + 1,1);
    
    for o = 1:length(Orders)
        %next states are computed as current inventory (o-1) + orders -
        %possible demands. 
        %
        % Note: no checking is performed for states > max_inv, so only pass
        % me valid orders (an non-negative demands)
        
        %-- setup
        total_units = o-1+Orders(o);
        
        %-- compute state transition matix
        next_s = total_units-demand;
        % Some next states will be negative indicating some unmet demand.
        % since we don't take backorders, we set these states to zero.
        next_s = max(next_s, 0);
        
        % now compute the transition probabilities by
        % first lumping the zero and negative inventories:
        prob = sum(demandprob(next_s <= 0));
        % then tack on the non-negative probabilities by
        %   first identifying the indexes for non-zero next states
        s_idx = next_s(next_s >0);
        %   then shifting to match the indexes starting at one for states
        s_idx = s_idx+1;
        % and finally matcing the indexes to the corresponding
        % probabilities
        prob(s_idx) = demandprob(next_s > 0);
        
        %and finally put this into the current row, and let the existing
        %zero padding take care of any impossible states
        P(o,1:length(prob)) = prob;

        % -- compute corresponding expected income
        %compute the sales & income for each demand
        sales = min(total_units, demand);
        
        %notice that when the next state is zero inventory, we must have
        %sold total_units, rather than demand. Hence we can use a similar
        %agregation via indexing that we used above
        incomes = income_fun(total_units);
        incomes(s_idx) = income_fun(sales(next_s > 0));
        
        ev_incomes(o) = sum(incomes .* prob);
    end
end % trans_matrix function