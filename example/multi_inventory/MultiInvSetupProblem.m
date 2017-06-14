function multi_inv_problem = MultiInvSetupProblem(multiinv_params, varargin)
% MULTIINVPROBLEMSETUP Setup dynamo problem structure for mult-inventory
% problem
%
% multi_inv_problem = MultiInvSetupProblem(multiinv_params)
%    Sets up a multi inventory problem structure for use with dynamo. This
%    includes assigning required and optional functions and other settings 
%
% multi_inv_problem = MultiInvSetupProblem(__, 'field_name', value, ...)
%    Selectively override defaults using arbitrary sets of key-value pairs.
%    Commonly used fields include:
%      'discount_rate'
%      'n_periods'
%
% multi_inv_problem = MultiInvProblemSetup(__, alt_problem_setup)
%    Selectively override defaults with values passed in as a structure or
%    key-value pair cell array
%
%
% Examples (and doctest):
%
% >> params = {'total_space', 3; 'unit_space', 1; 'p_demand', [0.25; 0.5; 0.25] };
% >> params = MultiInvSetupParams(params);
% >> multi_inv_problem = MultiInvSetupProblem(params)
%
% multi_inv_problem = 
% 
%   struct with fields:
% 
%                     params: [1×1 struct]
%              discount_rate: 0.1000
%                  n_periods: 3
%                  state_set: {1×4 cell}
%             fTerminalValue: @MultiInvTerminalValue
%               fDecisionSet: @MultiInvDecisionSet
%              fDecisionCost: @MultiInvDecisionCost
%             fDecisionApply: @MultiInvDecisionApply
%          decision_vfun_map: []
%               random_items: {[1×1 rpDiscreteSample]}
%                fRandomCost: @MultiInvRandomCost
%               fRandomApply: @MultiInvRandomApply
%           random_state_map: {[]}
%         fOpsBeforeDecision: []
%          fOpsAfterDecision: @MultiInvOps
%            fOpsAfterRandom: []
%     fCompareDecisionPolicy: []
%              fCompareValue: []
%             fMapState2Vfun: []
%           fOptimalDecision: []
%              fRandomSample: []
%               fRandomJoint: []
%
% >> multi_inv_problem = MultiInvSetupProblem(params, 'discount_rate', 0)
%
% multi_inv_problem = 
% 
%   struct with fields:
% 
%                     params: [1×1 struct]
%              discount_rate: 0
%                  n_periods: 3
%                  ****
%
% >> alt_as_cell = {'discount_rate', 0.2; 'n_periods', 5};
% >> multi_inv_problem = MultiInvSetupProblem(params, alt_as_cell)
%
% multi_inv_problem = 
% 
%   struct with fields:
% 
%                     params: [1×1 struct]
%              discount_rate: 0.2000
%                  n_periods: 5
%                  ****
%
% >> alt_as_struct.discount_rate = 0.03;
% >> alt_as_struct.n_periods = 4;
% >> multi_inv_problem = MultiInvSetupProblem(params, alt_as_struct)
%
% multi_inv_problem = 
% 
%   struct with fields:
% 
%                     params: [1×1 struct]
%              discount_rate: 0.0300
%                  n_periods: 4
%                  ****

% See also: MultiInvSetupParams

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   2  2017-06-14 06:33  BryanP      Doctests and fixed subfunctions 
%   1  2017-06-14 05:06  BryanP      Extracted from MultiInv_demo v18 


%===Now Setup general structure for MultiInv problem
multi_inv_problem_defaults = struct(...
        ...%Problem Setup
        'params',                   multiinv_params, ...   % Problem specific placeholder to be passed to all functions TODO add a struct before this
        'discount_rate',            0.1, ...
        'n_periods',                3, ...
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

multi_inv_problem = DefaultOpts(varargin, multi_inv_problem_defaults);
    
%===Create required ADP object instances
% NOTE: after problem structure defined because cell arrays get split
% within calls to struct()

%--Create (pre-decision) state space (state_set)
if isequal(multi_inv_problem.state_set, 'assign_later_to_avoid_error')
    multiinv_state_set = { setCombinWithLimits('',  false, multiinv_params.unit_space, multiinv_params.total_space)};
    % Since state is the same for all time periods we first create it and
    % then replicate it.
    % Notes: 
    %  -- b/c sets are handle classes, this only makes a shallow copy, which is 
    %     OK since we don't modify the sets during the algorithm runs.
    %  -- Need n_periods+1 to cover terminal value states
    multi_inv_problem.state_set = repmat( multiinv_state_set, 1, multi_inv_problem.n_periods+1);
end

%--Create random processes (random_items)
%
if isequal(multi_inv_problem.random_items, 'assign_later_to_avoid_error')
    multi_inv_problem.random_items = cell(1, multi_inv_problem.params.n_products);
    for p_idx = 1:multi_inv_problem.params.n_products
        vals = { ( 0:multi_inv_problem.params.max_inv(p_idx) )' };  % column vector within a cell array
        prob = multi_inv_problem.params.p_demand(p_idx);    %Already a cell array, just
        multi_inv_problem.random_items{p_idx} = rpDiscreteSample(vals, prob);
    end
end
if isequal(multi_inv_problem.random_state_map, 'assign_later_based_on_random_items')
    %Since DiscreteSample is independant of state, configure the
    %random_state_map with all empties
    multi_inv_problem.random_state_map = cell(size(multi_inv_problem.random_items));
end
