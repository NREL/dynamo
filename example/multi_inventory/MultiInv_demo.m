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
% Test using the small case and pre-computed optimal policy:
% >> small_opt_policy = {[0,0]		[2,4]	[2,4]	[1,2]; [0,1]		[2,3]	[2,3]	[0,0]; [0,2]		[2,2]	[2,2]	[0,0];  [0,3]		[2,1]	[0,0]	[0,0] ; [0,4]		[2,0]	[2,0]	[0,0] ; [0,5]		[2,0]	[2,0]	[0,0] ; [0,6]		[0,0]	[0,0]	[0,0];  [1,0]		[1,4]	[1,4]	[0,2] ; [1,1]		[1,3]	[1,3]	[0,0]; [1,2]		[0,0]	[0,0]	[0,0] }; 
% >> small_opt_value = 'TODO';
% >> [small_prob, small_result] = MultiInv_demo('small', 'dp');
%     Backward Induction DP
%         T=4 (terminal period): Done
%         T=3:Done: 44 states
%         T=2:Done: 44 states
%         T=1:Done: 44 states
%     Elapsed time is *** seconds.
%
% >> small_policy = cell2mat(small_result.dpbi_policy);
% >> small_opt_policy = cell2mat(small_opt_policy);
% >> isequal(small_opt_policy(:, 3:end), small_policy(1:10,:))
%    ans = 
%         1

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

    sbi_opt = { 'sbi_state_samples_per_time'            150    % Number of state samples per time period
                'sbi_decisions_per_sample'              20     % Number of decision samples per state
                'sbi_uncertain_samples_per_post'        10     % Number of random/uncertainty samples per time, used for all decisions
                
                'vfun_approx'                           'LocalRegr'
                };

end

% Build up additional derived fields. Also takes poisson p_demand and
% converts to a pdf in a cell array
multiinv_params = MultiInvParamsSetup(multiinv_params);
             
%===Now Setup general structure for MultiInv problem
multi_inv_problem = struct(...
        ...%Problem Setup
        'params',                   multiinv_params, ...   % Problem specific placeholder to be passed to all functions TODO add a struct before this
        'discount_rate',            multiinv_discount, ...
        'n_periods',                multiinv_n_periods, ...
        ...% State Related
        'state_set',                'assign_later_to_avoid_error', ...    % A cell vector of set objects capturing the pre-decision states. One per timestep TODO: support single item
        'fTerminalValue',           @MultiInvTerminalValue, ...      % Returns a terminal value for a list of terminal states
        ...% Decision Related
        'fDecisionSet',             @MultiInvDecisionSet, ...    % Returns a set object for a given pre-decision state
        'fDecisionCost',            @MultiInvDecisionCost, ...   % a function handle that returns the decision cost for a given list of decisions
        'fDecisionApply',           @MultiInvDecisionApply, ...  % Returns a post decision state list given a list of pre-decision states and decision
        'decision_vfun_map',        [], ...    % A map of indicies from decision dimensions to value function dimensions to allow adding decision cost and value function approximations
        ...% Uncertainty (aka Random) Related
        'random_items',             'assign_later_to_avoid_error', ... % A cell vector of RandProc objects,
        'fRandomCost',              @MultiInvRandomCost, ...  % Returns the cost associated with each random sample for the corresponding post-decision state
        'fRandomApply',             @MultiInvRandomApply, ... % Returns list of next pre-decision states given a post-decision state and list of random samples to apply
        'random_state_map',         'assign_later_based_on_random_items', ...    % A cell array map of indicies to extract state-tracked random process information (e.g. lattice). cell entry order matches random_items
        ...%Operations Cost Related: 0-3 may be defined as needed. If not defined, zero cost is assumed
        'fOpsBeforeDecision',       [], ... %Simulate operations and return operations costs based on pre-decision state (Optional) 
        'fOpsAfterDecision',        @MultiInvOps, ... % Simulate operations and return operations costs based on post-decision state and decision (Optional)
        'fOpsAfterRandom',          [], ... % Simulate operations and return operations costs based on post-decision state, decision, and uncertainty (before advancing time) (Optional)
        ...% Utility Functions (Optional).
        'fCompareDecisionPolicy',   [], ... % Compare two optimal policy decision sets (e.g. to check convergence)
        'fCompareValue',            [], ... % Compare two value function sets (e.g. to check convergence)
        'fMapState2Vfun',           [], ... % Map the full state space to an alternate value function space
        ...% Optional performance enhancements to replace general algorithms with custom designs.
        'fOptimalDecision',         [], ... % Find the least cost combination of Decision cost and Value Function for a given pre-decision state
        'fRandomSample',            [], ... % Replace the default sampling across all random_items RandProc objects
        'fRandomJoint',             [] ... % Replace the default joint distribution assembly across all random_items RandProc objects
        );

%===Create required ADP object instances
% NOTE: after problem structure defined because cell arrays get split
% within calls to struct()

%--Create (pre-decision) state space (state_set)
% TODO: replace with setCombinWithLimits
multiinv_state_set = { setCombinWithLimits('',  false, multiinv_params.unit_space, multiinv_params.total_space)};
% Since state is the same for all time periods we first create it and
% then replicate it.
% Notes: 
%  -- b/c sets are handle classes, this only makes a shallow copy, which is 
%     OK since we don't modify the sets during the algorithm runs.
%  -- Need n_periods+1 to cover terminal value states
multi_inv_problem.state_set = repmat( multiinv_state_set, 1, multiinv_n_periods+1);

%--Create random processes (random_items)
%
multi_inv_problem.random_items = cell(1, multi_inv_problem.params.n_products);
for p_idx = 1:multi_inv_problem.params.n_products
    vals = { ( 0:multi_inv_problem.params.max_inv(p_idx) )' };  % column vector within a cell array
    prob = multi_inv_problem.params.p_demand(p_idx);    %Already a cell array, just
    multi_inv_problem.random_items{p_idx} = rpDiscreteSample(vals, prob);
end
%Since DiscreteSample is independant of state, configure the
%random_state_map with all empties
multi_inv_problem.random_state_map = cell(size(multi_inv_problem.random_items));

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


% Old cruft below here. Could be starting point for other sized 

% %% ---- Single Item Treated as Multi Inventory
% % Common problem parameters
% max_space = 3;
% item_space = [1]; %#ok<NBRAK>
% n_periods = 4;
% p_demand = [0.1, 0.2, 0.4, 0.3]; % This is the probability for 0:n units of demand
% disc_rate = 0;
%
% % Pure single item DP for comparision
% fprintf('\nStarting simple single product inventory...')
% tic
% [simple_orders, simple_values] = inventory_dp(max_space, n_periods, p_demand, disc_rate)
% toc
%
% %% dummy "Multi-product" with only a single DP example "super simple single product case as multi-product"
% % Note: the results should match the results above... and they do!
%
% fprintf('\nStarting MultiInv framework for the same single product example...')
% simplemulti_params=MultiInvInit(max_space, item_space, {p_demand});
%
% tic
% [simplemulti_orders, simplemulti_values] = DP(n_periods, disc_rate, ...
%     @MultiInvState, @MultiInvDecision, @MultiInvTermValue, @MultiInvCost, ...
%     @MultiInvTransProb, simplemulti_params)
% toc
%
% %% Now run same single item problem as ADP
%
% % %For now comment out while waiting for 1-D function approx to be available
% %
% % %Setup basic problem structure
% % simple_adp = MultiInvInit(max_inv, item_space, {p_demand});
% % simple_adp.n_periods = n_periods;
% % simple_adp.disc_rate = disc_rate;
% % % Only needed for adpSBI simple_prob.dyn_var = rpDiscreteSample();
% % simple_adp.first_state = [0]; %#ok<NBRAK>
% %
% % Adapt other options for 2-D case below
%
% %% ---- Small Two Item Multi Inventory ---
% % Note this is a multi product inventory model with storage space
% % constraints and independent poisson demands.
%
% %% DP for Small 2-item
%
% max_space = 6;
% item_space = [2 3];
% p_demand = [1 2]; % Note, scalar per product treated as poisson lambda values
% n_periods = 3;
% disc_rate = 0.1;
%
% %Setup problem
% multi2_params=MultiInvInit(max_space, item_space, p_demand);
%
% fprintf('\nStarting small 2 product case...\n')
% tic
% [multi2_orders, multi2_values] = DP(n_periods, disc_rate, ...
%     @MultiInvState, @MultiInvDecision, @MultiInvTermValue, @MultiInvCost, ...
%     @MultiInvTransProb, multi2_params, 50)
% toc
%
% %% ADP (TD1) for small 2-item
%
% %Setup basic problem structure
% multi2_adp = MultiInvInit(max_space, item_space, p_demand);
% multi2_adp.n_periods = n_periods;
% multi2_adp.disc_rate = disc_rate;
% % Only needed for adpSBI simple_prob.dyn_var = rpDiscreteSample();
% multi2_adp.first_state = [0 0];
%
% % Define dimensions
% multi2_adp.dims.pre_state = length(item_space);
% multi2_adp.dims.post_state = length(item_space);
% multi2_adp.dims.vfun_state = length(item_space);
% multi2_adp.dims.decision = length(item_space);
%
% % Setup required functions
% multi2_adp.fSim = @MultiInvSim;
% multi2_adp.fApplyDscn = @MultiInvApplyDecision;
% multi2_adp.fDecision = @MultiInvDecisionList;
% multi2_adp.fFirstState = @MultiInvFirstState;
%
% %Setup state randomizers for bootstrap
% multi2_adp.rand.state = rpDiscreteSample({MultiInvState('all', [], multi2_adp)});
% multi2_adp.rand.decision = @MultiInvDecisionSample;
%
% % Now setup ADP options
% % WARNING: The current version of faInterp does not effectively support
% % machine learning needs.
% adp_opts = {
%                 'vfun_approx'           'Interp'
%                 'boot_iter_per_t'       5
%                 'boot_sample_state'     true
%                 'boot_enable_full_sim'  false
%                 'max_iter'              500
%                 'explore_iter'          100
%                 'bkps_use_updated_vfun' false
%                 'bkps_abort'            true
%                 'fix_rand'              true
%            };
%
%
% % Finally Run ADP
% tic
% multi2_adp_results = adpTD1(multi2_adp, adp_opts)
% toc
%
% %% ---- Medium: still small enough to solve with DP ---
% % Note this is a multi product inventory model with storage space
% % constraints and independent poisson demands.
%
% %% DP for Medium 3-item
%
% max_space = 50;
% item_space = [1 3 5];
% p_demand = [30 10 5]; % Note, scalar per product treated as poisson lambda values
% n_periods = 3;
% disc_rate = 0.1;
%
% %Setup problem
% fprintf('\nSetting up medium 3 product case...\n')
% tic
% med3_params=MultiInvInit(max_space, item_space, p_demand);
% toc
%
% fprintf('\nStarting medium 3 product case...\n')
% tic
% [med3_orders, med3_values] = DP(n_periods, disc_rate, ...
%     @MultiInvState, @MultiInvDecision, @MultiInvTermValue, @MultiInvCost, ...
%     @MultiInvTransProb, med3_params, 50)
% toc
%
% %% ADP (TD1) for medium 3-item
%
% %Setup basic problem structure (building off of previous setup which can be
% %very slow for such large problems
% med3_adp = med3_params;
% med3_adp.n_periods = n_periods;
% med3_adp.disc_rate = disc_rate;
% % Only needed for adpSBI simple_prob.dyn_var = rpDiscreteSample();
%
% % Define dimensions
% med3_adp.dims.pre_state = length(item_space);
% med3_adp.dims.post_state = length(item_space);
% med3_adp.dims.vfun_state = length(item_space);
% med3_adp.dims.decision = length(item_space);
% med3_adp.first_state = zeros(1, med3_adp.dims.pre_state) ;
%
% % Setup required functions
% med3_adp.fSim = @MultiInvSim;
% med3_adp.fApplyDscn = @MultiInvApplyDecision;
% med3_adp.fDecision = @MultiInvDecisionList;
% med3_adp.fFirstState = @MultiInvFirstState;
%
% %Setup state randomizers for bootstrap
% med3_adp.rand.state = rpDiscreteSample({MultiInvState('all', [], med3_adp)});
% med3_adp.rand.decision = @MultiInvDecisionSample;
%
% % Now setup ADP options
% % WARNING: The current version of faInterp does not effectively support
% % machine learning needs.
% adp_opts = {
%                 'vfun_approx'           'Interp'
%                 'boot_iter_per_t'       50
%                 'boot_sample_state'     true
%                 'boot_enable_full_sim'  false
%                 'max_iter'              1e10    %Let time limit be binding
%                 'max_time_sec'          300
%                 'explore_iter'          100
%                 'bkps_use_updated_vfun' false
%                 'bkps_abort'            true
%                 'fix_rand'              true
%            };
%
%
% med3_adp_results = adpTD1(med3_adp, adp_opts)
