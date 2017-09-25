function [orders, results, value_function_approx] = ...
            inventory_adp(N, n_iter, inv_p, adp, restart_soln, ref_policy, tol)
%INVENTORY_ADP One product stochastic inventory approx. dynamic prog. (ADP)
%
% Usage: [orders, results,, value_function_approx] = ...
%            inventory_adp(N, n_iter, inv_p, adp, restart_soln, ref_policy, tol)
%
% Description: An approximate dynamic program (ADP) of the stochastic inventory 
% problem, as outlined in Puterman Section 3.2.2.  The program solves for 
% the optimal unit ordering under uncertain customer demand.  
%
% IMPORTANT: the value function approximation has changed. Details below
%
% ADP Details:
%   - Value function approximation: CAVE, piecewise linear based on slopes 
%   - Sampling: Monte Carlo pure-explore for bootstrap (default 100 trials)
%     then pure-exploitation (take apparent optimum) for later iterations
%   - Learning/Update: double-pass 
% 
% Inputs [with defaults in brackets];
%  N          Optional (skip with []) number of decision epochs or time
%              periods. This setting overrides that from inv_prob. Included
%              for backward compatability.
%  n_iter     Number of ADP iterations. If tol is also specified, this is
%              treated as a maximum number of iterations. If a restart
%              solution is provided, only n_iter - restart_soln.n additional
%              iterations, if any, are run [required]
%  inv_prob   Structure defining the inventory problem parameters possible
%              fields include [defaults in brackets]:
%                max_inv         maximum inventory that can be held [50]
%                n_periods       number of decision epochs or time periods [10]
%                dr              discount rate [5% = 0.05] 
%                lambda          poisson demand parameter [20]
%                init_state      initial inventory state [0]
%              plus special fields for cost function components:
%                order_cost_ea   Per unit ordering cost [2/ea]
%                order_cost_fix  Fixed cost per order [4/order]
%                order_min_qty   Minimum order quantity [1]
%                sales_price     Income on sales [8/ea]
%                hold_cost       Cost to keep in inventory [1/ea]
%                term_reward     Value of any inventory in final period [0/ea]
%
%  adp_setup  Structure describing the ADP approaches to take in solving
%              the problem. Valid fields [defaults in brackets]:
%               vfun_approx    Value function approximation, options:
%                    -'CAVE'     Piecewise linear concave [default]
%               adp.cave_intercept
%                              For CAVE only: boolean value to capture the 
%                               actual value function magnitude (true) or
%                               just the slopes [false, default]
%               bootstrap      Number of initial samples to monte carlo
%                               using our value approximation [100]
%               stepsize       Algorithm for determining stochastic
%                               gradient stepsize in learning. Options:
%                    -'1overN'       Yep, 1/n
%                    -'constant'     1 for first sample, constant specified
%                                    by step_opt for all others
%                    -'harmonic'     Generalized harmonic stepsize [default]
%                    - a string (e.g. 'Method') corresponding to a
%                      function of the form:
%                           step = ssMethod(n, old_step, step_opt)
%               step_opt       Optional, step size algorithm parameter.
%                               Can be a scalar, vector, structure, etc.
%                               depending on chosen stepsize function.
%                               [default = 10, corresponding to the 'a'
%                               harmonic parameter]
%               fix_rand       Use a repeatable stream of random numbers.
%                               Equivalant to starting with the same seed
%                               but with the fancy MATLAB rand [false]
%               plot           Option to plot the first specified number of 
%                               periods of value function [5]
%               trace          Store iteration trace time & tol data [false] 
%               inv4err        Inventory states to use for err calcs [0:10]            
%
%  restart_soln An existing adp solution from which to resume. Format is
%                the same as results output [default: build new solution]
%  ref_policy   Reference optimal order policy (such as from inventory_dp) to
%                compare for trace & converence tol. Can either be a full
%                array with one col per period and 1 row per state or a
%                vector for each period with the optimal postdecision state
%                [default: empty]
%  tol          Convergence tolerance for stopping. Treated as a relative
%                convergence unless ref_soln provided in which case it
%                provides an absolute tolerance. Additional convergence
%                parameters can be specified in the adp_setup struct.
%                [default: run until n_iter]
%
% Outputs:
%  orders     array of orders in the optimal policy. One column per time
%              period and one row per pre-decision inventory state
%  results    struct array containing solution results include the 
%              post-decision state value function approximation:
%                n          number of iterations
%                step       final step size (for iterative stepsizes)
%                vfun       value function approximation (1 row per period)
%                   for CAVE:
%                       vfun.states     states (x_values)
%                       vfun.slopes     corresponding slopes
%                postdec    optimal post decision state (ignoring order cost)
%                trace      iteration tracing, captures time, convergence,
%                            and other metrics about the solution process.
%                            Formated as a struct array with one row per
%                            iteration
% value_function_approx 
%             cell array containing the actual states and values for the 
%              value function approximation. There are N rows in the array
%              with columns as follow: 
%                  col 1: states
%                  col 2: value function magnitude for each state
%
% 5 key elements for this problem (Powell Chapter 5):
% State Variables:
% Decision Variables:
% Exogenous Information Processes:
% Transition Function:
% Objective Function:
%
%IMPORTANT: the value function approximations have changed. They are now
% the true post decision state and hence may not appear to have a finite
% maximum. This is correct according to the ADP forms of Bellman's equation
% found in Powell (2007). When ordering costs are introduced, the curve is
% effectively tilted such that there is a finite maximum as seen in the
% orders matrix
% ----------------------------------------------------------------------- %

% Implementation Notes:
%  - We shorten inv_prob to inv_p and adp_setup to adp to keep the code compact
%
% HISTORY
% ver     date    time        who      changes made
% ---  ---------- -----  ------------- ---------------------------------------
% 1-7  2010-06           MortW, NidhiS Initial debugging iterations
%                        & JenM
%   8  2010-08-07 11:41  BryanP        Removed version from filenames
%   9  2010-09-21 16:38  BryanP        Added option to ignore/use intercept
%  10  2010-09-21 16:47  JenM          Made boostrap # a parameter
%  11  2010-11-01 22:04  BryanP        Modularized inputs & defaults
%  12  2010-11-03 03:45  BryanP        MAJOR overhaul to match Powell:
%                                        - remove contribution from postdec value function
%                                        - actually optimize orders at each step
%                                        - explicit optimal policy computation
%                                        - Misc code streamlining & comments
%  13  2010-11-03 23:00  BryanP        Implemented more efficient find_opt_order
%  14  2010-11-04 01:30  BryanP        Corrected bug in find_opt_order opt_net_value 
%  15  2010-11-04 01:30  BryanP        Corrected treatment of discounting 
%                                        & location of stepsize
%  16  2010-11-05 04:15  BryanP        Clean up of extraneous code and updated
%                                      documentation for new output & options
%  17  2010-11-05 14:15  BryanP        MAJOR Enhancements:
%                                       - Ability to trace time & error at each iteration
%                                       - Ability to resume a solution
%                                       - Converted results output to struct
%  18  2010-11-07 22:15  BryanP        Refinements:
%                                       - Run until converged to ref solution
%                                       - Only include specified states in error calcs 
%                                       - Compute estimated time without trace calculations
%                                       - Store both relative and absolute errors 
%  19  2012-02-29 15:55  BryanP        Refinements:
%                                       - converge_check_n 
%                                       - Additional display: setup & not converged 
%                                       - String based stepsize functions

%===========================%
%      Handle Defaults      %
%===========================%

% If no problem parameters are specified, create blank structures to allow
% later code to set the defaults
if nargin < 3
    inv_p = struct();
end
if nargin < 4
    adp = struct();
end

%Define a cell array of defaults so we can simply loop through rather than
%having a long chain of if statements
inv_prob_defaults = ...
    {%field     %default
     'max_inv'    50        % maximum inventory that can be held
     'n_periods'  10        % number of decision epochs or time periods
     'dr'         0.05      % discount rate
	 'lambda'     20        % poisson demand parameter
     'init_state' 0         % initial inventory state
	 'order_cost_ea'    2
	 'order_cost_fix'   4
     'order_min_qty'    1
     'sales_price'      8
	 'hold_cost'        1
	 'term_reward'      0   % per unit value of any inventory in final period
    };

%If the structure field is not set or is blank, use the default
for idx = 1:length(inv_prob_defaults)
    fieldname = inv_prob_defaults{idx, 1};
    if not(isfield(inv_p, fieldname)) || isempty(inv_p.(fieldname))
        inv_p.(fieldname) = inv_prob_defaults{idx, 2};
    end
end
        
adp_setup_defaults = ...
    {%field         %default
	 'vfun_approx'    'CAVE'      %Value function approximation
     'cave_smooth'    [5, 5]      %CAVE: initial smoothing interval
     'cave_intercept' false       %CAVE: Capture actual value not just slopes 
     'ds'             1           %delta state for numeric derivative
	 'bootstrap'      100         %Number of monte carlo bootstrap samples
	 'stepsize'       'harmonic'  %Stocastic gradient stepsise algorithm 
	 'step_opt'       10          %step size algorithm parameter(s)
	 'fix_rand'       false       %(don't) use a repeatable stream of random numbers
     'plot'           5           %plot the first specified number of periods of value function
     'trace'          false       %Store iteration trace time & tol data
     'inv4err'        0:10        %Inventory states to use for err calcs
     'converge_check_n' 10        %Interval for convergence checks 
    };

for idx = 1:length(adp_setup_defaults)
    fieldname = adp_setup_defaults{idx, 1};
    if not(isfield(adp, fieldname)) || isempty(adp.(fieldname))
        adp.(fieldname) = adp_setup_defaults{idx, 2};
    end
end

%Force the value of N for number of periods
if not(isempty(N))
    inv_p.n_periods = N;
end

% Force trace to be on if iterating until converged to tolerence so that we
% compute error each pass
if nargin > 6
    adp.trace = true;
end

% Decide which type of solution error to report based on presence or absense
% of a ref_policy
if nargin >= 6 && not(isempty(ref_policy))
    find_abs_err = true;
    is_converged=false;
else
    find_abs_err = false;
end
 
%===========================%
%      Initialization       %
%===========================%

% Compute Additional parameters
limits = [0 inv_p.max_inv]; % limits for function approximation (used by CAVE)
inv_p.disc_factor = 1 - inv_p.dr; % discount factor, put in inv_p for helper function

% code to get consistent random numbers during debugging
if adp.fix_rand
    disp('inventory_adp using fixed random number stream')
	rand_num_sequence = RandStream('mt19937ar');
	RandStream.setDefaultStream(rand_num_sequence)
end

%======== Initialize output storage =======
%Use restart solution if provided
if nargin >5 && not(isempty(restart_soln))
    n_start = restart_soln.n;
    step = restart_soln.step;
    
    %Allocate extra trace storage
    if adp.trace
    	results.trace.time(n_start:n_iter) = zeros(n_iter-n_start, 1);
    	results.trace.rel_err(n_start:n_iter) = zeros(n_iter-n_start, 1);
    end

else
    %Otherwise intialize a new problem
    n_start = 1;
    step = [];              % initialize current step size to override the Simulink function step()
    
    %And allocate space for the solution (each is a column struct vector)
    results.vfun.states = num2cell(zeros (inv_p.n_periods+1, 1));
    results.vfun.slopes = num2cell(zeros (inv_p.n_periods+1, 1));
    if adp.trace
        results.trace.time = NaN * ones(n_iter, 1);
        results.trace.time_wo_trace = NaN * ones(n_iter, 1);
        results.trace.rel_err = NaN * ones(n_iter, 1);
        if find_abs_err
        	results.trace.abs_err = NaN * ones(n_iter, 1);
        end
    end
end

% initialize data storage vectors for each iteration
predec_state = zeros(inv_p.n_periods+1,1);
order = zeros(inv_p.n_periods+1,1);
postdec_state = zeros(inv_p.n_periods+1,1);

% initalize trace
% Note that using t= tic allows for running our own nested timer without
% affecting an outer timer for the entire function
% TODO use last trace value as starting point
if adp.trace
    disp('inventory_adp trace enabled, run times will be much slower')
    extra_trace_time = 0;
    trace_t = tic;
end

%Display setup
fprintf('ADP Inventory: v_fun=CAVE, step=%s, conv_check=%d, max_iter=%d\n', ...
    adp.stepsize, find_abs_err, n_iter)
if not(any(strcmp(adp.stepsize, {'1overN', 'constant', 'harmonic'})))
    adp.stepsize = str2func(['ss' adp.stepsize]);
end

%===========================%
%   Run ADP (double pass)   %
%===========================%

% ADP outer iteration loop
for n = n_start:n_iter
    
    %gradually reduce smoothing interval
    %TODO only for CAVE
    if (n < n_iter/10)
        smooth = adp.cave_smooth;
    elseif (n < n_iter/2)
        smooth = adp.cave_smooth/2;
    else
        smooth = adp.cave_smooth/5;
    end
     
    % randomly sample poisson demands for entire time (but don't let the
    % algorithm have any foresight)
    demand = poissrnd(inv_p.lambda, inv_p.n_periods+1,1);
    
    % initialize existing inventory aka first period (inv_p.n_periods=1) 
    %pre-decision state variable
    predec_state(1) = inv_p.init_state;
   
    
    %===== FORWARD PASS =======
    % loop forward over decision periods
    for t = 1:inv_p.n_periods
                
        % randomly choose order quantity for bootstrap iterations,
        % then optimize to find the best decision for remaining iterations
        if (n <= adp.bootstrap)
            max_order = inv_p.max_inv - predec_state(t);

            %order(t) = unidrnd(max_order + 1) - 1;
            %if (predec_state(t) >= inv_p.max_inv)
            %    order(t) = max_order;
            %end

            order(t) = poissrnd(inv_p.lambda);
            if (order(t) > max_order)
                order(t) = max_order;
            end
        else
            % TODO: use the references directly in the interpolation
            % Start by extracting the current approximation for the next
            % time period
            next_vfun_states = results.vfun.states{t+1};
            next_vfun_slopes = results.vfun.slopes{t+1};     

            % Now call our helper function to find the optimal order
            % TODO: handle non-CAVE
            order(t)  = find_opt_order(predec_state(t), next_vfun_states, next_vfun_slopes, inv_p);
            
        end
        
        % compute postdec_state for this order
        postdec_state(t) = predec_state(t) + order(t);

        % Handle random sales and compute next period predec_state
        predec_state(t+1) = max((postdec_state(t) - demand(t)), 0);

        %TODO: store value function approximation? to save multiple calls to
        %FunFromSlope

        %TODO: and/or alternatively go ahead and compute the slopes for CAVE?
    
    end
    
    %Handle (deterministic) stepsizes
    % Do this here rather than in the backward pass to avoid duplicate
    % computations and to work properly with iterative stepsize
    % computation such as McClain's
    %
    % Note: changed 'alpha' to 'step' to aviod conflict with
    % plotting transparency option function alpha()
    if isa(adp.stepsize, 'function_handle')
        step = adp.stepsize(n, step, adp.step_opt);
    else
        if strcmp(adp.stepsize, 'harmonic')
            step = adp.step_opt / (adp.step_opt + n - 1);
        elseif strcmp(adp.stepsize, 'constant')
            if n == 1
                step = 1;
            else
                step = adp.step_opt;
            end
        elseif strcmp(adp.stepsize, '1overN')
            step = 1/n;
        else
            error('DP_examples:inventory:inventory_adp:stepsize', ...
                    ['Unknown stepsize function, ' adp.stepsize]);
        end
    end
  
    %===== BACKWARD PASS =======
    % Work backward to update the value function
    for t = inv_p.n_periods:-1:1
        
        %TODO: directly reference these results?
        % pull out current time period cell array vectors (again) to use
        states = results.vfun.states{t};
        slopes = results.vfun.slopes{t};     
        
        if t == inv_p.n_periods
            % Note: the final period (n_period + 1) value function is known
            % and is a simple linear function with slope equal terminal
            % reward so it is easy to represent in a CAVE like way
            next_states = 0;
            next_slopes = inv_p.term_reward; 
            %next_intercept = 0;
        else
            next_states = results.vfun.states{t+1};
            next_slopes = results.vfun.slopes{t+1}; 
        end
     
        %Find the actual post decision value function and associated slopes
        postdec_for_slopes = min(inv_p.max_inv, max(0, postdec_state(t) + [-adp.ds 0 adp.ds]));

        %Compute realized costs for these orders/post decision states
        for s_idx = 1:length(postdec_for_slopes);
            % Compute actual sample value as sales reward - cost of resulting 
            % optimal order + corresponding postdec_state value for next period
            [opt_order, next_net_value] = ...
                find_opt_order(max(0,postdec_for_slopes(s_idx) - demand(t)), next_states, next_slopes, inv_p); %#ok<ASGLU>

            % Note: We need to adjust the income from sales to make up for
            % the fact that we will discount the entire value function
            % approximation, including the sales income component, when
            % finding the optimal order in the forward pass
            sample_values(s_idx) = inv_p.sales_price/inv_p.disc_factor * min(postdec_for_slopes(s_idx),demand(t)) ...
                                        + inv_p.disc_factor * next_net_value;
        end
                
        %TODO: CAVE specific
        % compute 1 X 2 vector of sample slopes
        sample_slopes = diff(sample_values)./adp.ds;
            
        %%% sample_slopes

        % update value function approximation with current sample

        %ToDo: make CAVE Specific
        [new_states new_slopes] = CaveUpdateStep(postdec_state(t), sample_slopes, states, slopes, limits, step, smooth);

        
        % replace existing cell array vectors with new values for states and slopes
        results.vfun.states{t} = new_states;
        results.vfun.slopes{t} = new_slopes;
    end
    
    %Handle trace
    if adp.trace
        results.trace.time(n) = toc(trace_t);
        results.trace.time_wo_trace(n) = results.trace.time(n) - extra_trace_time;
        orders = find_optimal_policy(results, inv_p, adp.inv4err);
        
        if n == 1
            old_orders = zeros(size(orders));
        end
        cur_policy_err = orders - old_orders;
        results.trace.rel_err(n) = norm(cur_policy_err, 'fro')/norm(orders, 'fro');
        old_orders = orders;
        
        if find_abs_err && mod(n,adp.converge_check_n)==0
            cur_policy_err = orders - ref_policy(adp.inv4err+1,:);
            results.trace.abs_err(n) = norm(cur_policy_err, 'fro')/norm(ref_policy, 'fro');
            
            %If a convergence tolerance is defined
            if nargin > 6
                if results.trace.abs_err(n) < tol
                    fprintf('Yay, policy CONVERGED to %g (<%g) in only %d iterations\n',results.trace.abs_err(n), tol, n)
                    is_converged=true;
                    break
                end
            end
        end
        %Compute the cumlative extra time for "w/o trace time"
        extra_trace_time = extra_trace_time + toc(trace_t) - results.trace.time(n);
    end
end

if find_abs_err && not(is_converged)
    fprintf('DID NOT CONVERGE after %d iterations\n', n)
end

%===========================%
%      Format Output        %
%===========================%

orders = find_optimal_policy(results, inv_p);

results.n = n;
results.step = step;

%===========================%
%      Handle plotting      %
%===========================%

%If needed for output or for plotting, compute the value function
%approximation magnitudes.
if nargout > 2 || isfield(adp,'plot') && adp.plot > 0
	value_function_approx = cell(inv_p.n_periods,2);
    
    % calculate value function approximation for each time period and store
    for t = 1:inv_p.n_periods

        %ToDo: make CAVE Specific
        [total_values states] = FunFromSlope(results.vfun.states{t}, results.vfun.slopes{t}, inv_p.max_inv);
        value_function_approx{t,2} = total_values;
        value_function_approx{t,1} = states;


    end
end

% Do the plotting as requred
if isfield(adp,'plot') && adp.plot > 0
	linespecs = {'-.r' '-.g' '-.b' '-.c' '-.m' '-.y' '-.k', 'r' 'g' 'b' 'c' 'm' 'y' 'k','.r' '.g' '.b' '.c' '.m' '.y' '.k' };

    hold on
    for t = 1:adp.plot

        plot(value_function_approx{t,1},value_function_approx{t,2}, linespecs{t});

    end
    title(sprintf('first %d periods of the value function approximation', adp.plot))
    xlabel('inventory')
    legend (int2str((1:adp.plot)'));
end

end %main function

%===========================%
%      Helper Functions     %
%===========================%

%TODO: Vectorize?
function [opt_order, opt_net_value] = ...
        find_opt_order(predec_state, next_vfun_states, next_vfun_slopes, inv_p)

    % Helper function to find optimal order for a piecewise linear value
    % function approximation

    % We have an affine (linear + offset) opt_order cost, which could result in
    % a non-convex optimization for the Bellman equation directly. So
    % instead we break into two parts: the best non-zero order and the no
    % order case and compare the results
                        
    % No order: only have to pay holding costs
    %TODO check this vs next value function approximation
    no_order_value = - inv_p.hold_cost * predec_state ...
                        + inv_p.disc_factor * InterpFromSlope(next_vfun_states, next_vfun_slopes, predec_state);
    
    % if we are already at the maximum inventory, we can't order any more
    % so just return this value
    if predec_state == inv_p.max_inv
        opt_order = 0;
        opt_net_value = no_order_value;
        return
    end

	% For the opt_order case, we take advantage of the fact that our per unit
	% opt_order costs can easily be added to the value function within slope
	% space. This is b/c the piece-wise linear slope of the bellman
	% function will simply be the sum of the marginal decision cost (per
	% unit hold + opt_order costs) + the piecewise linear future value slopes
    %
    % We also have to:
    %   - offset the value by the opt_order's fixed cost
    %   - limit our search space to the range [predec_state, max_inv]
    %
    % Implementation details:
    %   - We add 1 to the predec_state (and later to the opt_order) since
    %     the minimum opt_order is 1 and we treat the no opt_order case separately.
    %     For positive fixed opt_order cost this is actually not neccessary,
    %     but this implementation is robust to funky ordering incentives,
    %     etc.
    
    % Adjust slopes to reflect ordering costs
    %
    % Note1: the slopes here are only valid for post-decision states greater
    % than our current pre-decision state. This is because the order cost
    % does not apply to items currently in inventory (don't have to order
    % the pre-decision state inventory). Computing the correct set of
    % slopes would typically involve inserting a new vertex (at the
    % predecision state) into the approximation, which seems messy.
    % Instead, we:
    %   - Take care to only use the portion of the approximation where the
    %     shape (slopes) are OK. This is the portion corresponding to valid
    %     post-decision states higher than our predec_state. This
    %     restriction is handled by checking for non-positive orders, and
    %   - Adjust the entire approximation up by the ammount which we
    %     over-estimated the ordering costs (in_stock_value_adjust) to
    %     provide correct value magnitudes
    %
    % Note2: For discounting, we only want to apply discounting to the
    % expected future value from next_vfun. Here we take advantage of the
    % fact that a * y for y = m*x+b --> a*m*x + b*x. This allows us to
    % discount the future value function simply by discounting the slope
    % values Within this space the intercept, b, is ignored so discounting
    % does not affect it.
    net_value_slopes = inv_p.disc_factor * next_vfun_slopes - inv_p.order_cost_ea - inv_p.hold_cost;
    in_stock_value_adjust = inv_p.order_cost_ea * predec_state;
    
    % Find optimal postdecision state
    %
    % Note: to get accurate value estimates we shift the y_values by
    % in_stock_adjustment as described above
    [opt_postdec_state, opt_order_value] = ...
        MaxFromSlope(next_vfun_states, net_value_slopes, inv_p.max_inv, ...
                        in_stock_value_adjust);
    %Force state to integer
    opt_postdec_state = floor(opt_postdec_state);
    
	% And find the corresponding optimal order
    opt_order = opt_postdec_state - predec_state;
	
    % If reaching the optimal postdecision state requires a non-positive
    % order, we can't get there, so use the minimum order state instead
    %
    % Notes:
    %   - to get accurate value estimates we shift the y_values by
    %     in_stock_adjustment as described above.
    %   - We have to compute the adjusted value function for the actual
    %     postdec_state this order will take us to.
    if opt_order < inv_p.order_min_qty
        opt_order = inv_p.order_min_qty;
        opt_order_value = InterpFromSlope(next_vfun_states, net_value_slopes,...
                                            predec_state + opt_order, ...
                                            in_stock_value_adjust);
    end
	opt_order_value = opt_order_value - inv_p.order_cost_fix;
    
    %Pick best case:
    if no_order_value >= opt_order_value
        % Best not to order
        opt_order = 0;
        opt_net_value = no_order_value;
    else
        % Order away
        opt_net_value = opt_order_value;
    end
        
end

%compute the optimal policy (orders) for specified (pre-decision) states and time periods
function orders = find_optimal_policy(results, inv_p, inv_states)
    if nargin <3
        inv_states = 0:inv_p.max_inv;
    end
    orders = zeros(length(inv_states), inv_p.n_periods);

    for t = 1:inv_p.n_periods
        for s = inv_states
            orders(s+1, t) = find_opt_order(s, results.vfun.states{t}, results.vfun.slopes{t}, inv_p); 
        end
    end
end
