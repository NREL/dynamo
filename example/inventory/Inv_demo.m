% Demonstration for:
% modularized Dynamic Programming implementation of SINGLE PRODUCT Inventory
% Originally Developed as part of MIT ESD.862 Final Class Project
% Bryan Palmintier
% Spring 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   5  2017-09-24 22:20  BryanP      Overhaul complete for new problem structure
%   4  2017-09-23 23:40  BryanP      Renamed Inv_demo and adopt ideas from MultiInv_demo v19 and MultiInvSetupProblem v3 
%   3  2010-08-07 11:50  BryanP      Removed MultiInv to focus only on Inv
%   2  2010-05-26 23:50  BryanP      MultiInvCombin tests
% many 2010              BryanP      original code, additions not tracked

%% ========== Question 1 - Inventory Problem (Finite Horizon) ==========

%% --- Part 1: Repeat Putterman 3.2.2
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

fprintf('\n\n----- Running small inventory problem (Putterman 3.2.2) -----\n')

p_demand = [0.25, 0.5, 0.25]';
max_inv = 3;
%get problem defaults & compute poisson demand distribution
simple_inv_params = InvSetupParams(max_inv,p_demand);

simple_inv_problem = struct(...
        ...%Problem Setup
        'params',                   simple_inv_params, ...   % Problem specific placeholder to be passed to all functions TODO add a struct before this
        'discount_rate',            0, ...
        'n_periods',                3, ...
        ...% State Related
        'state_set',                {{setList( (0:max_inv)' )}}, ...
        'fTerminalValue',           @InvTerminalValue, ...      % Returns a terminal value for a list of terminal states
        ...% Decision Related
        'fDecisionSet',             @InvDecisionSet, ...    % Returns a set object for a given pre-decision state
        'fDecisionCost',            @InvDecisionCost, ...   % a function handle that returns the decision cost for a given list of decisions
        'fDecisionApply',           @InvDecisionApply, ...  % Returns a post decision state list given a list of pre-decision states and decision
        'decision_vfun_map',        [], ...    % A map of indicies from decision dimensions to value function dimensions to allow adding decision cost and value function approximations
        ...% Uncertainty (aka Random) Related
        'random_items',             {{rpDiscreteSample( {(0:2)'}, {p_demand} ) }}, ... % A cell vector of RandProc objects,
        'fRandomCost',              @InvRandomCost, ...  % Returns the cost associated with each random sample for the corresponding post-decision state
        'fRandomApply',             @InvRandomApply, ... % Returns list of next pre-decision states given a post-decision state and list of random samples to apply
        'random_state_map',         {{[]}}, ...    % A cell array map of indicies to extract state-tracked random process information (e.g. lattice). cell entry order matches random_items. Empty if no state tracking
        ...%Operations Cost Related: 0-3 may be defined as needed. If not defined, zero cost is assumed
        'fOpsBeforeDecision',       [], ... %Simulate operations and return operations costs based on pre-decision state (Optional) 
        'fOpsAfterDecision',        @InvOps, ... % Simulate operations and return operations costs based on post-decision state and decision (Optional)
        'fOpsAfterRandom',          [] ... % Simulate operations and return operations costs based on post-decision state, decision, and uncertainty (before advancing time) (Optional)
    );

tic
[SimpleDP_Results] = dpBI(simple_inv_problem);
toc

sbi_opt = { 'sbi_state_samples_per_time'            20    % Number of state samples per time period
            'sbi_decisions_per_sample'              5    % Number of decision samples per state
            'sbi_uncertain_samples_per_post'        5     % Number of random/uncertainty samples per time, used for all decisions

            'vfun_approx'                           'LocalAvg'
            };

tic
[SimpleADPsbi_Results] = adpSBI(simple_inv_problem, sbi_opt);
toc


%% --- Part 2: Larger system with various Poisson demands
time_steps = 15;
max_lambda = 20;
max_inv = 10;

% Adapt common portions of problem from small example
poisson_inv_problem = simple_inv_problem;
poisson_inv_problem.n_periods = time_steps;
poisson_inv_problem.state_set = {setList( (0:max_inv)' )};

fprintf('\nRunning modularized DP for inventory problem with %d different lambda values\n', max_lambda)
clear PoissonDP_Results
tic
for lambda = 1:max_lambda
    %revise problem defaults & compute poisson demand distribution
    this_poisson_inv_params = InvSetupParams(max_inv, lambda);
    
    poisson_inv_problem.params = this_poisson_inv_params;
    poisson_inv_problem.random_items = ...
        {rpDiscreteSample( {(0:max_inv)'}, {this_poisson_inv_params.p_demand} ) };


    [PoissonDP_Results{lambda}] = dpBI(poisson_inv_problem); %#ok<SAGROW>
end
fprintf('Total Poisson sensitivity dpBI ')
toc

%% --- Part 3: Add Discounting
time_steps = 15;
max_inv = 10;
lambda = 5;

disc_rates = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5];

% Adapt common portions of problem from latest poisson example
inv_disc_problem = simple_inv_problem;
inv_disc_problem.n_periods = time_steps;
inv_disc_problem.state_set = {setList( (0:max_inv)' )};
inv_disc_problem.params = InvSetupParams(max_inv, lambda);
inv_disc_problem.random_items = ...
    {rpDiscreteSample( {(0:max_inv)'}, {inv_disc_problem.params.p_demand} ) };


clear DiscDP_Results
tic
for dr = 1:length(disc_rates)
    inv_disc_problem.discount_rate = disc_rates(dr);
	DiscDP_Results{dr} = dpBI(poisson_inv_problem); %#ok<SAGROW>
end
fprintf('Total discount rate sensitivity dpBI ')
toc


