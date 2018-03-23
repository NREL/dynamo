function simple_storage_pv_problem = SimpleStoragePvSetupProblem(simple_storage_pv_params, varargin)
% SIMPLESTORAGEPVPROBLEMSETUP Setup dynamo problem structure for simple storage + PV problem 
%
% simple_storage_pv_problem = SimpleStoragePvSetupProblem(simple_storage_pv_params)
%    Sets up a storage-PV problem structure for use with dynamo. This
%    includes assigning required and optional functions and other settings 
%
% simple_storage_pv_problem = SimpleStoragePvSetupProblem(__, 'field_name', value, ...)
%    Selectively override defaults using arbitrary sets of key-value pairs.
%    Commonly used fields include:
%      'discount_rate'
%
% simple_storage_pv_problem = SimpleStoragePvProblemSetup(__, alt_problem_setup)
%    Selectively override defaults with values passed in as a structure or
%    key-value pair cell array
%
%
% Examples (and doctest):
%

% See also: SimpleStoragePvSetupParams

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   2  2017-07-16 21:24  BryanP      Completed initial parameter, etc. definition 
%   1  2017-07-12 12:43  BryanP      MultiInvSetupProblem v2


%===Now Setup general structure for SimpleStoragePv problem
simple_storage_pv_problem_defaults = struct(...
        ...%Problem Setup
        'params',                   simple_storage_pv_params, ...   % Problem specific placeholder to be passed to all functions TODO add a struct before this
        'discount_rate',            0.04, ... 
        'n_periods',                simple_storage_pv_params.n_periods, ... % Number of decision periods
        ...% State Related
        'state_set',                'assign_later_to_avoid_error', ...    % A cell vector of set objects capturing the pre-decision states. One per timestep TODO: support single item
        'fTerminalValue',           @SimpleStoragePvTerminalValue, ...      % Returns a terminal value for a list of terminal states
        ...% Decision Related
        'fDecisionSet',             @SimpleStoragePvDecisionSet, ...    % Returns a set object for a given pre-decision state
        'fDecisionCost',            @SimpleStoragePvDecisionValue, ...  % a function handle that returns the decision cost for a given list of decisions
        'fDecisionApply',           @SimpleStoragePvDecisionApply, ...  % Returns a post decision state list given a list of pre-decision states and decision
        'decision_vfun_map',        [], ...    % A map of indicies from decision dimensions to value function dimensions to allow adding decision cost and value function approximations
        ...% Uncertainty (aka Random) Related
        'random_items',             'assign_later_to_avoid_error', ... % A cell vector of RandProc objects,
        'fRandomCost',              @SimpleStoragePvRandomValue, ...  % Returns the cost associated with each random sample for the corresponding post-decision state
        'fRandomApply',             @SimpleStoragePvRandomApply, ... % Returns list of next pre-decision states given a post-decision state and list of random samples to apply
        'random_state_map',         'assign_later_based_on_random_items', ...    % A cell array map of indicies to extract state-tracked random process information (e.g. lattice). cell entry order matches random_items
        ...%Operations Cost Related: 0-3 may be defined as needed. If not defined, zero cost is assumed
        'fOpsBeforeDecision',       @SimpleStoragePvOps, ... %Simulate operations and return operations costs based on pre-decision state (Optional) 
        'fOpsAfterDecision',        @SimpleStoragePvOps, ... % Simulate operations and return operations costs based on post-decision state and decision (Optional)
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

simple_storage_pv_problem = DefaultOpts(varargin, simple_storage_pv_problem_defaults);
    
%===Create required ADP object instances
% NOTE: after problem structure defined because cell arrays get split
% within calls to struct()

% Extract parameters for easier access
params = simple_storage_pv_problem.params;

%--Create (pre-decision) state space (state_set)
% Order: pv_pen storage_kW storage_kWh
if isequal(simple_storage_pv_problem.state_set, 'assign_later_to_avoid_error')
    simple_storage_pv_problem.state_set = cell(1, params.n_periods);
    simple_storage_pv_states = SimpleStoragePvSetupStates(simple_storage_pv_problem.params);
    for t = 1:params.n_periods + 1 %Include terminal period
        simple_storage_pv_problem.state_set{t} = setList(simple_storage_pv_states{t});
    end
end

%--Create random processes (random_items)
%
if isequal(simple_storage_pv_problem.random_items, 'assign_later_to_avoid_error')
    simple_storage_pv_problem.random_items = cell(1);
    simple_storage_pv_problem.random_items{1} = rpTransMatrix(params.pv_state_set, params.pv_trans_set);
end
if isequal(simple_storage_pv_problem.random_state_map, 'assign_later_based_on_random_items')
    %ID which column in state correspond to the current pv penetration
    % Note: random_state_map is a cell vector listing the sets of indicies
    % that correspond to matching random_items entries
    simple_storage_pv_problem.random_state_map = { 1 };
end
