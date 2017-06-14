function [multi_inv_problem, results] = MultiInv_demo(varargin)
% MULTIINV_DEMO Demo script showing simple use of adp clases of algorithms
%
% multi_inv_problem = MultiInv_demo()
%   set up full example multi_inventory problem structure
%
% [multi_inv_problem, multiinv_params] = MultiInv_demo()
%   also separately return inner parameters
%
% MultiInv_demo(scenario_str)
%   specify which scenario to setup. Valid options are:
%     'single'   small, single inventory problem (Putterman 3.2.3)
%     'small'    A 2 dimensional multi-inventory that runs in a few seconds for DP 
%     'medium'   A 3 dimensional multi-inventory that runs in about a minute for DP 
%
% MultiInv_demo(__, 'dp')
% MultiInv_demo(__, 'dpBI')
%   Run traditional backward induction on example MultiInv problem
%
% MultiInv_demo(__, 'sbi')
% MultiInv_demo(__, 'adpSBI')
%   Run sampled backward induction ADP algorithm on example MultiInv problem
%
%
% %% Test DP algorithm using the small case and compare to pre-computed optimal policy:
% >> small_opt_policy_cell = {[0,0]		[2,4]	[2,4]	[1,2]; [0,1]		[2,3]	[2,3]	[0,0]; [0,2]		[2,2]	[2,2]	[0,0];  [0,3]		[2,1]	[0,0]	[0,0] ; [0,4]		[2,0]	[2,0]	[0,0] ; [0,5]		[2,0]	[2,0]	[0,0] ; [0,6]		[0,0]	[0,0]	[0,0];  [1,0]		[1,4]	[1,4]	[0,2] ; [1,1]		[1,3]	[1,3]	[0,0]; [1,2]		[0,0]	[0,0]	[0,0] }; 
%
% >> small_opt_value = 'TODO';
% >> [small_prob, small_dp_result] = MultiInv_demo('small', 'dp');
%     Backward Induction DP
%         T=4 (terminal period): Done
%         T=3:Done: 44 states
%         T=2:Done: 44 states
%         T=1:Done: 44 states
%     Elapsed time is *** seconds.
%
% >> small_dp_policy = cell2mat(small_dp_result.dpbi_policy);
% >> small_opt_policy = cell2mat(small_opt_policy_cell);
% >> isequal(small_opt_policy(:, 3:end), small_dp_policy(1:10,:))
%    ans = 
%         logical
%         1
%
%
% %% And test adpSBI algorithm using small case:
% >> rng('default')
% >> [~, small_sbi_result] = MultiInv_demo('small', 'sbi');
% Sampled Backward Induction ([0.1] ksamples/period)
%     Creating empty post-decision value functions (LocalRegr)
%     T=4 (terminal period): Done
%     T=3:S........................................100
%     T=2:S........................................100
%     T=1:S........................................100
% Warning: Multiple initial states defined, using first in list 
% > In adpSBI (line ***)
%   **** 
% Elapsed time is *** seconds.
%
% >> small_sbi_result
% 
% small_sbi_result = 
% 
%   struct with fields:
% 
%     first_decision: [2 4]
%          objective: 16.1959
%          post_vfun: [1×4 faLocalRegr]
%            adp_opt: [1×1 struct]
% %
% >> isequal(small_sbi_result.first_decision, small_opt_policy_cell{1,2})
% 
% ans =
%      logical
%      1


% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  18  2017-04-26 17:32  BryanP      Added many pre-sized examples  
%  17  2017-04-26 05:24  BryanP      BUGFIX: use setCombinWithLimits (finally available)  
%  16  2017-04-09 22:42  BryanP      Added fRandomJoint for DP support 
%  15  2017-04-03 10:28  BryanP      Convert to a function to allow selection of algorithms to try 
%  14  2016-12-01 22:15  BryanP      Use [] rather than NaN for unused functions 
%  13  2016-11-10 13:15  BryanP      Changed terminal values to mixed non-zero values 
%  12  2016-11-10 12:45  BryanP      create state_set entry for n_periods+1 to cover terminal value states
%  11  2016-11-10 12:25  BryanP      Added blocks to run DP & ADP 
%  10  2016-11-10 11:45  BryanP      Renamed write out "Terminal" in fTerminalValue 
%   9  2016-10-21 16:31  BryanP      Add terminal value function 
%   8  2016-10-21 15:31  BryanP      FIX stepsize in set definition 
%   7  2016-10-07 01:31  BryanP      Problem definition done? Added random_items 
%   6  2016-10-06 11:31  BryanP      Build up params and sets 
%   5  2016-09-16 14:35  BryanP      Complete draft problem structure
%   4  2016-09-01 11:15  BryanP      Rework for new problem structure
%   3  2016-05-01 05:23  BryanP      ADP runs (poor results). And added huge case
%   2  2016-04    11:23  BryanP      Added basic 2 product case in DP
%   1  2016-04-14 11:23  BryanP      Adapted from MultiInvDP_scratchpad v5

%% Setup problem specific values
if any(strcmpi('single', varargin))
    % Reproduce Putterman 3.2.2
    %
    % Create "single" test case if requested. Which runs a single product
    % example as a degenerate multi-product example.
    % With old MultiInvDP_scratchpad this ran in ~5.2sec on MacBook Pro mid-2014
    %
    % Solution
    % Optimal Policy:
    % state  t= 1   2   3
    %   0		3	2	0
    %   1		0	0	0
    %   2		0	0	0
    %   3		0	0	0
    %
    % Optimal Values:
    %  state    t= 1  	  2     3   4 
    %    0      4.1875   2      0	0
    %    1      8.0625	 6.25	5	0
    %    2      12.125	 10     6	0
    %    3      14.1875	 10.5	5	0
    multiinv_n_periods = 3;
    multiinv_discount = 0;
    multiinv_params = { %REQUIRED: Warehouse and product space, size, & demand configuration
                        'total_space'       3      % Max space in warehouse
                        'unit_space'        1       % Space per item
                        'p_demand'          [0.25; 0.5; 0.25] % Probability vector (pmf) of sales per quantity per item. If scalar per item: treated as poisson demand, cell with columns otherwise
                       % OPTIONAL: Pricing assumptions (can use defaults).
                       % Note: costs are negative
                        'order_cost'        -4      % Flat cost to place any order (in additon to per unit costs)
                        'term_unit_val'     0       % Value of product at end of forecast horizon
                        'unit_cost'         -2      % Ordering cost per unit. If scalar, assumes same cost for all products
                        'hold_cost'         -1      % Cost to keep in warehouse per unit per time. If scalar, assumes same cost for all products
                        'sales_price'       8       % Sales price per unit. If scalar, assumes same cost for all products
                     };    
                 
	sbi_opt = { 'sbi_state_samples_per_time'            20    % Number of state samples per time period
                'sbi_decisions_per_sample'              5    % Number of decision samples per state
                'sbi_uncertain_samples_per_post'        5     % Number of random/uncertainty samples per time, used for all decisions
                
                'vfun_approx'                           'LocalAvg'
                };

elseif any(strcmpi('medium', varargin)) || any(strcmpi('med', varargin))
    % Create "medium" test case if requested. Matches the old
    % MultiInvDP_scratchpad "medium" case which runs in 17.4sec on MacBook Pro
    % mid-2014 model
    %
    % Optimal Solution (top 10 entries)
    % >> optimal_policy = {
    % %state      t= 1       2       3       4       5
    % [0,0,0]     [2,3,5]	[2,3,5]	[2,3,5]	[2,3,5]	[1,2,3]
    % [0,1,0]     [2,2,5]	[2,2,5]	[2,2,5]	[2,2,5]	[1,1,3]
    % [0,2,0]     [2,1,5]	[2,1,5]	[2,1,5]	[2,1,5]	[1,0,3]
    % [0,3,0]     [2,0,5]	[2,0,5]	[2,0,5]	[2,0,5]	[1,0,3]
    % [0,4,0]     [2,0,4]	[2,0,4]	[2,0,4]	[2,0,4]	[1,0,3]
    % [0,5,0]     [1,0,3]	[1,0,3]	[1,0,3]	[1,0,3]	[1,0,3]
    % [0,6,0]     [0,0,2]	[0,0,2]	[0,0,2]	[0,0,2]	[0,0,2]
    % [1,0,0]     [1,3,5]	[1,3,5]	[1,3,5]	[1,3,5]	[0,2,3]
    % [1,1,0]     [1,2,5]	[1,2,5]	[1,2,5]	[1,2,5]	[0,1,3]
    % [1,2,0]     [1,1,5]	[1,1,5]	[1,1,5]	[1,1,5]	[0,0,3] };
    %
    % >> optimal_value = {
    % % state    t= 1  	  2       3       4 	  5      6
    % [0,0,0]     73.1	60.4	46.4	30.8	13.3	0.0
    % [0,1,0]     75.1	62.4	48.4	32.8	15.3	0.0
    % [0,2,0]     77.1	64.4	50.4	34.8	17.3	0.0
    % [0,3,0]     79.1	66.4	52.4	36.8	18.9	0.0
    % [0,4,0]     80.7	68.1	54.0	38.2	19.1	0.0
    % [0,5,0]     79.8	67.1	53.0	36.8	18.5	0.0
    % [0,6,0]     73.7	61.0	46.8	30.5	14.0	0.0
    % [1,0,0]     75.1	62.4	48.4	32.8	15.3	0.0
    % [1,1,0]     77.1	64.4	50.4	34.8	17.3	0.0
    % [1,2,0]     79.1	66.4	52.4	36.8	19.3	0.0 };
    %
    % >> [med_prob, med_result] = MultiInv_demo('medium', 'dp');
    % ***
    % >> med_policy = cell2mat(med_result);
    % >> optimal_policy = cell2mat(optimal_policy);
    % >> isequal(optimal_policy(:, 4:end), med_policy(1:10,:))
    %
    % ans =
    %      1

    multiinv_n_periods = 5;
    multiinv_discount = 0.1;
    multiinv_params = { %REQUIRED: Warehouse and product space, size, & demand configuration
                        'total_space'       20      % Max space in warehouse
                        'unit_space'        [2 3 1] % Space per item
                        'p_demand'          [1 2 3] % Probability vector (pmf) of sales per quantity per item. If scalar per item: treated as poisson demand, cell with columns otherwise
                       % OPTIONAL: Pricing assumptions (can use defaults).
                       % Note: costs are negative
                        'order_cost'        -4      % Flat cost to place any order (in additon to per unit costs)
                        'term_unit_val'     0       % Value of product at end of forecast horizon
                        'unit_cost'         -2      % Ordering cost per unit. If scalar, assumes same cost for all products
                        'hold_cost'         -1      % Cost to keep in warehouse per unit per time. If scalar, assumes same cost for all products
                        'sales_price'       8       % Sales price per unit. If scalar, assumes same cost for all products
                     };    

	sbi_opt = { 'sbi_state_samples_per_time'            200    % Number of state samples per time period
                'sbi_decisions_per_sample'              30    % Number of decision samples per state
                'sbi_uncertain_samples_per_post'        30     % Number of random/uncertainty samples per time, used for all decisions
                
                'vfun_approx'                           'LocalAvg'
                };

elseif any(strcmpi('large', varargin))
    % NOTE: Not well tested
    
    multiinv_n_periods = 3;
    multiinv_discount = 0.1;
    multiinv_params = { %REQUIRED: Warehouse and product space, size, & demand configuration
                        'total_space'       50      % Max space in warehouse
                        'unit_space'        [1 3 5] % Space per item
                        'p_demand'          [30 10 5] % Probability vector (pmf) of sales per quantity per tiem. If scalar per item: treated as poisson demand, cell with columns otherwise
                       % OPTIONAL: Pricing assumptions (can use defaults).
                       % Note: costs are negative
                        'order_cost'        -4      % Flat cost to place any order (in additon to per unit costs)
                        'term_unit_val'     [1 2 3] % Value of product at end of forecast horizon
                        'unit_cost'         -2      % Ordering cost per unit. If scalar, assumes same cost for all products
                        'hold_cost'         -1      % Cost to keep in warehouse per unit per time. If scalar, assumes same cost for all products
                        'sales_price'       8       % Sales price per unit. If scalar, assumes same cost for all products
                     };    

	sbi_opt = { 'sbi_state_samples_per_time'            400    % Number of state samples per time period
                'sbi_decisions_per_sample'              80    % Number of decision samples per state
                'sbi_uncertain_samples_per_post'        80     % Number of random/uncertainty samples per time, used for all decisions
                
                'vfun_approx'                           'LocalAvg'
                };
else
    % DEFAULT TO SMALL PROBLEM
    %
    % Create "small" test case if requested. Also run as old
    % MultiInvDP_scratchpad "small" case in ~0.5sec on MacBook Pro mid-2014
    % model
    %
    % Note: does not run when using doctest
    %
    % Optimal Solution (top 10 entries)
    % >> small_opt_policy = {
    % % state    t= 1       2       3
    % [0,0]		[2,4]	[2,4]	[1,2]
    % [0,1]		[2,3]	[2,3]	[0,0]
    % [0,2]		[2,2]	[2,2]	[0,0]
    % [0,3]		[2,1]	[0,0]	[0,0]
    % [0,4]		[2,0]	[2,0]	[0,0]
    % [0,5]		[2,0]	[2,0]	[0,0]
    % [0,6]		[0,0]	[0,0]	[0,0]
    % [1,0]		[1,4]	[1,4]	[0,2]
    % [1,1]		[1,3]	[1,3]	[0,0]
    % [1,2]		[0,0]	[0,0]	[0,0] };
    %
    % >> small_opt_policy = 'TODO';
    % >> [small_prob, small_result] = MultiInv_demo('small', 'dp');
    % ...
    % >> small_policy = cell2mat(small_result.dpbi_policy);
    % >> small_opt_policy = cell2mat(small_opt_policy);
    % >> isequal(small_opt_policy(:, 3:end), small_policy(1:10,:))

    multiinv_n_periods = 3;
    multiinv_discount = 0.1;
    multiinv_params = { %REQUIRED: Warehouse and product space, size, & demand configuration
                        'total_space'       20      % Max space in warehouse
                        'unit_space'        [2 3] % Space per item
                        'p_demand'          [1 2] % Probability vector (pmf) of sales per quantity per item. If scalar per item: treated as poisson demand, cell with columns otherwise
                       % OPTIONAL: Pricing assumptions (can use defaults).
                       % Note: costs are negative
                        'order_cost'        -4      % Flat cost to place any order (in additon to per unit costs)
                        'term_unit_val'     0       % Value of product at end of forecast horizon
                        'unit_cost'         -2      % Ordering cost per unit. If scalar, assumes same cost for all products
                        'hold_cost'         -1      % Cost to keep in warehouse per unit per time. If scalar, assumes same cost for all products
                        'sales_price'       8       % Sales price per unit. If scalar, assumes same cost for all products
                     };    

    sbi_opt = { 'sbi_state_samples_per_time'           100    % Number of state samples per time period
                'sbi_decisions_per_sample'              20     % Number of decision samples per state
                'sbi_uncertain_samples_per_post'        10     % Number of random/uncertainty samples per time, used for all decisions
                
                'vfun_approx'                           'LocalRegr'
                };

end

% Build up additional derived fields. Also takes poisson p_demand and
% converts to a pdf in a cell array
multiinv_params = MultiInvSetupParams(multiinv_params);
multi_inv_problem = MultiInvSetupProblem(multiinv_params, ...
    'disc_rate', multiinv_discount, 'n_periods', multiinv_n_periods);
             
%% Now that setup is complete, let's run the specified examples
if nargin < 1
    return
end


%% Run DP (Backward Induction) algorithm
if any(strcmpi('dp', varargin)) || any(strcmpi('dpBI', varargin))
    multi_inv_problem_dp = multi_inv_problem;

    tic
    results = dpBI(multi_inv_problem_dp);
    toc
    %TODO: display summary
    return
end

%% Run ADP Sample Backward Induction algorithm
if any(strcmpi('sbi', varargin)) || any(strcmpi('adpSBI', varargin))
    multi_inv_problem_sbi = multi_inv_problem;

    tic
    results = adpSBI(multi_inv_problem_sbi, sbi_opt);
    toc
    %TODO: display summary
    return
end

%% Run ADP Double Pass (Temporal Difference, Lambda=1) algorithm
if any(strcmpi('td1', varargin)) || any(strcmpi('adpTD1', varargin))
    multi_inv_problem_td1 = multi_inv_problem;
    tic
    results = adpTD1(multi_inv_problem_td1);
    toc
    %TODO: display summary
    return
end

