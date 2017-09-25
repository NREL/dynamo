function [results value_function_approx all_sample_total_values] = inventory_adp_regress(N, iter)

% Authors:
% Nidhi R. Santen (MIT ESD) and Maria Teresa Pena Alcaraz (IIT, Comillas
% University), with help from Mort Webster (MIT ESD)
% October 2010
%
% Description: An approximate dynamic program of the stochastic inventory 
% problem, as outlined in Puterman Section 3.2.2.  The program solves for 
% the optimal unit ordering under uncertain customer demand.  A recursive least-squares regression method 
% is used to approximate the value function.  A double-pass
% algorithm is used.
%
% N is the number of decision epochs (here, that is also time periods)
% iter is the number of samples
%
% results is a cell array, containing the regression coefficients and optimal
% state for each time period.  
% There are N+1 rows in the array.  
% regression coefficients are in column 1
% optimal state is in column 2
%
% value_function_approx is a cell array, containing the states and values
% for the value function approximation.
% There are N rows in the array (state values and value function values for
% each time period)
%
%
% 5 key elements for this problem (Powell Chapter 5):
% State Variables:
% Decision Variables:
% Exogenous Information Processes:
% Transition Function:
% Objective Function:
%
% ----------------------------------------------------------------------- %
 
% % code to get consistent random numbers during debugging
% r = RandStream('mt19937ar');
% RandStream.setDefaultStream(r)
 
% vector with the characteristics of the plots 
linespecs = {'-.r' '-.g' '-.b' '-.c' '-.m' '-.y' '-.k' '-.r' '-.g' '-.b' ...
            '-.c' '-.m' '-.y' '-.k' '-.r' '-.g' '-.b' '-.c' '-.m' '-.y' ...
            '-.k' '-.r' '-.g' '-.b' '-.c' '-.m' '-.y' '-.k' '-.r' '-.g' ...
            '-.b' '-.c' '-.m' '-.y' '-.k'};
 
% initialize parameters
MAX = 15;         % MAX is the maximum inventory that can be held
dr = 0.05;        % discount rate
delta = 1 - dr;   % discount factor
lambda = 5;       % poisson demand parameter 
bootstrap = 500;  % number of samples for bootstrap
 
% initialize vector of all possible states for value function
% approximiation
all_states = (0:MAX);
 
% initialize storage cell array for regression coefficients and optimal_state
results = cell(N+1,2);  
 
% initialize storage cell array for value function approximation
value_function_approx = cell(N,2);
    
% initialize bootstrap sample_total_value storage cell array
bootstrap_sample_total_values = cell(N,2);     
 
% store all sample_total_values for debugging/polynomial fitting
all_sample_total_values = cell(N,2);     
 
%initialize storage cell array to track coeffiient updates over iterations
theta_updates = cell(iter-bootstrap,N);
 
% initialize sample reward vector (for i = 1 only)
sample_reward = zeros(N+1,1);
        
% initialize cell array to store most recent B for each time period in the recursive regression
B_holding = cell(N,1);
        
% initialize data storage vectors for forward pass for each iteration
predec_state = zeros(N+1,1);
order = zeros(N+1,1);
postdec_state = zeros(N+1,1);
 
% samples through states
for i = 1:iter
     
    % randomly sample poisson demands   
    demand = poissrnd(lambda, N+1,1);
    % demand = 2 * ones(N+1,1);
    
    % initialize first period (N=1) pre-decision state variable
    predec_state(1) = 0;
    
    % loop FORWARD over decision periods
    for t = 1:N+1       
        % randomly choose order quantity for first 3 samples (bootstrap),
        % then use order decision that matches the optimal (postdecision) state for remaining samples
        if (i <= bootstrap)                             % this bootstrap needs to match the bootstraps below b/c an optimal state does 
                                                        % not exist until after the first regression
            max_order = MAX - predec_state(t);
            order(t) = unidrnd(max_order + 1) - 1;
            if (predec_state(t) >= MAX)
                order(t) = 0;
            end
        else
            opt_state = results{t,2};                   % lookup opt_state from stored value in results array
            if (predec_state(t) < opt_state)
                order(t) = opt_state - predec_state(t);
            else
                order(t) = 0;
            end
        end
 
        % compute postdec_state for this order
        postdec_state(t) = predec_state(t) + order(t);
               
        % compute next period predec_state
        if (t < N+1)
            predec_state(t+1) = max((postdec_state(t) - demand(t)), 0);
        end 
    end
 
    % loop BACKWARDS over decision periods
    for t = N:-1:1  
              
        % first, compute current period's reward
        sample_reward(t) = inventory_cost_regress(postdec_state(t), order(t), demand(t)); 
        
        % calculate next period's value
        if (i <= bootstrap)
            if (t == N)
                sample_next_total_value = 0;
            else        
                sample_next_total_value = sample_total_value;
            end
        else
            if (t == N)
                sample_next_total_value = 0;
            else
                sample_next_total_value = interp1q(value_function_approx{t+1,1}', value_function_approx{t+1,2}', postdec_state(t+1));
            end
        end
                        
        % calculate Bellman Value, save
        sample_total_value = sample_reward(t) + (delta * sample_next_total_value);
        
        if (i <= bootstrap)
            bootstrap_sample_total_values{t,2}(i) = sample_total_value;
            bootstrap_sample_total_values{t,1}(i) = postdec_state(t);
        end
        
        all_sample_total_values{t,2}(i) = sample_total_value;
        all_sample_total_values{t,1}(i) = postdec_state(t);
 
        % now, run least-squares regression, save vector of 3 coefficients,
        % and calculate current value function approximation
        if (i >= bootstrap)
            if (i == bootstrap)
                % pull out needed parameters to use for this time period
                xsamp = all_sample_total_values{t,1};      
                ysamp = all_sample_total_values{t,2};
 
                % define new parameters for regression
                xsamp2 = xsamp .^ 2;
                xsamp3 = xsamp .^ 3;
                xsamp4 = xsamp .^ 4;
 
                % Now, run a regular OLS regression
                xinsub = [xsamp4' xsamp3' xsamp2' xsamp' ones(i,1)];
                ysampsub = ysamp';
                theta = regress(ysampsub, xinsub);                        % regression on full bootstrap
                theta_updates{i,t} = theta;                               % keeping track of starting coefficients                 
                results{t,1} = theta;
                
            else
                % now, start recursive least squares for the rest of the samples
                theta = results{t,1};                % pull out the most recent coefficients to update
 
                if (i == bootstrap + 1)              
                    B = inv(xinsub' * xinsub);       % for each time period, take the inverse of B one time (the original B)
                else
                    B = B_holding{t,1};
                end
                x_new = [(postdec_state(t))^4 (postdec_state(t))^3 (postdec_state(t))^2 postdec_state(t) 1]';     
                y_hat = x_new' * theta;          % this is the prediction of the new y based on the current approximation of the coefficients
                y_new = sample_total_value;      % this is the y value of the new sample
                epsilon = y_hat - y_new;         % this is the residual (current prediction - observation)
                gamma = 1 + x_new' * B * x_new;
                H = (1/gamma) * B;
                theta_new = theta - (H * x_new * epsilon);      
                theta_updates{i,t} = theta_new;       % keeping track of coefficients over iterations          
                results{t,1} = theta_new;             % here, the old regression coefficients is replaced with the new coefficients
                B_new = B - (1/gamma) * (B * x_new * x_new' * B);
                B_holding{t,1} = B_new;
            end
 
            % calculate new value function approximation
            new_total_values = ((results{t,1}(1) .* (all_states .^ 4)) + (results{t,1}(2) .* (all_states .^ 3)) + (results{t,1}(3) .* (all_states .^ 2)) + (results{t,1}(4) .* (all_states)) + results{t,1}(5));        
                              
            [C,I] = max(new_total_values);
            opt_state = all_states(I(1));         
            results{t,2} = opt_state;
            
            % save current value function approximation 
            value_function_approx{t,1} = all_states;   
            value_function_approx{t,2} = new_total_values;
        end
    end
end
 
for t = 1:N
     % Plot all the samples we have
     plot(all_sample_total_values{t,1},all_sample_total_values{t,2},'.');
     hold on
     % Plot all the value of the approximate function
     plot(value_function_approx{t,1},value_function_approx{t,2}, linespecs{t});    
end
 
end
