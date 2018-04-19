function HEMS_problem = HEMS_setup_problem( HEMS_params,varargin )
%function HEMS_problem = HEMS_setup_problem( HEMS_params,varargin )
%Adapted from MultiInv by BryanP
%Author: Li Wang
HEMS_problem_defaults = struct(...
        ...%Problem Setup
        'initial_state',            [78],...
        'params',                   HEMS_params, ...   % Problem specific placeholder to be passed to all functions TODO add a struct before this
        'discount_rate',            0, ...
        'n_periods',                24, ...
        ...% State Related
        'state_set',                'assign_later_to_avoid_error', ...    % A cell vector of set objects capturing the pre-decision states. One per timestep TODO: support single item
        'fTerminalValue',           @HEMS_TerminalValue, ...      % Returns a terminal value for a list of terminal states
        ...% Decision Related
        'fDecisionSet',             @HEMS_DecisionSet, ...    % Returns a set object for a given pre-decision state
        'fDecisionCost',            @HEMS_DecisionCost, ...   % a function handle that returns the decision cost for a given list of decisions
        'fDecisionApply',           @HEMS_DecisionApply, ...  % Returns a post decision state list given a list of pre-decision states and decision
        'decision_vfun_map',        [], ...    % A map of indicies from decision dimensions to value function dimensions to allow adding decision cost and value function approximations
        ...% Uncertainty (aka Random) Related
        'random_items',             'assign_later_to_avoid_error', ... % A cell vector of RandProc objects,
        'fRandomCost',              @HEMS_RandomCost, ...  % Returns the cost associated with each random sample for the corresponding post-decision state
        'fRandomApply',             @HEMS_RandomApply, ... % Returns list of next pre-decision states given a post-decision state and list of random samples to apply
        'random_state_map',         [], ...    % A cell array map of indicies to extract state-tracked random process information (e.g. lattice). cell entry order matches random_items
        ...%Operations Cost Related: 0-3 may be defined as needed. If not defined, zero cost is assumed
        'fOpsBeforeDecision',       [], ... %Simulate operations and return operations costs based on pre-decision state (Optional) 
        'fOpsAfterDecision',        [], ... % Simulate operations and return operations costs based on post-decision state and decision (Optional)
        'fOpsAfterRandom',          @HEMS_Ops, ... % calculate the discomfort cost
        ...% Utility Functions (Optional).
        'fCompareDecisionPolicy',   [], ... % Compare two optimal policy decision sets (e.g. to check convergence)
        'fCompareValue',            [], ... % Compare two value function sets (e.g. to check convergence)
        'fMapState2Vfun',           [], ... % Map the full state space to an alternate value function space
        ...% Optional performance enhancements to replace general algorithms with custom designs.
        'fOptimalDecision',         [], ... % Find the least cost combination of Decision cost and Value Function for a given pre-decision state
        'fRandomSample',            [], ... % Replace the default sampling across all random_items RandProc objects
        'fRandomJoint',             [] ... % Replace the default joint distribution assembly across all random_items RandProc objects
        );
HEMS_problem = DefaultOpts(varargin, HEMS_problem_defaults);

%--Create (pre-decision) state space (state_set)
if isequal(HEMS_problem.state_set, 'assign_later_to_avoid_error')
        HEMS_problem.state_set = cell(1,HEMS_problem.n_periods+1);
        cap_mat = cell2mat(HEMS_params.appliance_range);
        min = cap_mat(:,1)'./HEMS_params.unit_space;
        max = cap_mat(:,2)'./HEMS_params.unit_space;
       
        HEMS_problem.state_set{1} = setList(HEMS_problem.initial_state);
        for i = 2: HEMS_problem.n_periods+1
            HEMS_problem.state_set{i} = setCombinWithLimits('',  true, HEMS_params.unit_space, HEMS_params.total_space,max,0,min);
        end
        
        %assume pre-decision state spaces are same for all time period.
        %HEMS_state_set = { setCombinWithLimits('',  false, HEMS_params.unit_space, HEMS_params.total_space,max,0,min)};
        %HEMS_problem.state_set = repmat( HEMS_state_set, 1, HEMS_problem.n_periods+1); 
    
end
% Create uncertainty object for each appliance in every timeslots.
if isequal(HEMS_problem.random_items, 'assign_later_to_avoid_error')
     HEMS_problem.random_items = cell(1, 1);
     timeslots = HEMS_problem.params.time_slots;
     vals = cell(1,timeslots);
     prob = cell(1,timeslots);
     for p_idx = 1:1
         for t_idx = 1:timeslots
            s_range_mat = cell2mat(HEMS_params.stochastic_range);
            s_min = s_range_mat(p_idx,1);
            s_max = s_range_mat(p_idx,2);
            s_unit = HEMS_params.unit_space(p_idx);
            vals{t_idx} = (s_min:s_unit:s_max)';
            prob{t_idx} = HEMS_params.stochastic_prob{t_idx,p_idx};
         end
         HEMS_problem.random_items{p_idx} = rpDiscreteSample(vals, prob);
     end
end

