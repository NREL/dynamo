function results = adpSBI(problem, adp_opt, post_vfun)
% adpSBI Sampled Backward Induction for estimated ADP and intialization
%
% results = adpSBI(problem, adp_opt, post_vfun)
%
% Notes: 
% -- This method is currently only applicable to problem classes where it is
% possible to sample the post decision state without requiring a full,
% simulation derived state (ie it does not work for partially hidden
% Markov/quasi-Markov processes).
% -- post_vfun and fn inputs primarily provided for use 
% 
%
% Required problem attributes
%   n_periods
%   disc_rate
%   dyn_var list
%
% Required problem methods:
%   preToPost: takes decision list and pre_dec state and returns a list
%     of post decision states. This also  allows the user to map to a
%     reduced dimensionality
%   contrib(state_list, t): not counting decision costs
%
% Optional problem methods/functions
%   pre-dec state sample function
%   multi-contribution: compute multiple contribution functions at once.
%     this also allows user to work with a contribution function
%     approximation
%
% Required value function methods
%   update value function: takes a vector of (value based) states and their
%     corresponding value
%
% Adp options:
%   Sample function... pure Monte Carlo, hypercube, etc.
%
% Required sub-functions:

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  26  2017-04-03 11:12  BryanP      Working with new problem structure 
%  25  2017-04-02 06:00  BryanP      Made unique state collapse for ops costs optional and clean-up operations cost handling 
%  24  2017-03-06 06:00  BryanP      Added after decision ops computation and updated psuedocode  
%  23  2016-12-xx        BryanP      Further updates, nearly working with new problem structure 
%  22  2016-07-06 16:55  BryanP      Partial rework for new problem structure 
%  21  2012-07-06 16:55  BryanP      Added time to PostToVfun calls
%  20  2012-07-03 14:55  BryanP      Renamed to adpSBI (was adpSampBackInd)
%  19  2012-06-20/24     BryanP      Convert (par)for loop "Core" functions from nested to true sub-functions (easier debug)
%  18  2012-06-20 10:15  BryanP      Improved handling of bad (too sparse) lookup including sbi_min_ok_for_ev
%  17  2012-06-08 12:55  BryanP      New adp defaults: decfun_params
%  16  2012-05-14 12:45  BryanP      Added adp prefix to name
%  15  2012-05-11 16:25  BryanP      Fix objective to include discounting and decision contribution
%  14  2012-05-08 14:05  BryanP      Renamed from PreToPost to ApplyDscn
%  13  2012-05-07 11:22  BryanP      Add dimension names to vfun
%  12  2012-05-05 23:55  BryanP      Add adp to results
%  11  2012-04-22        BryanP      BUGFIX: future costs. AND param/return for nested functions
%  10  2012-04-21 17:20  BryanP      Rely on faLocalRegr for auto-expand of neighborhood
%   9  2012-04-19 21:45  BryanP      Separate postdec from vfun state spaces
%   8  2012-04-16 17:25  BryanP      Use generic optimal decision function
%   7  2012-04-07 11:25  BryanP      Added verbose
%   6  2012-04           BryanP      Problem specific sample + Refine CapPlan decisions
%   5  2012-03-26 16:45  BryanP      Use FuncApprox family
%   4  2012-03-20 16:00  BryanP      Use SampleNdRange for default sample objects
%   3  2012-03-19        BryanP      Hack in ncw expansion problem
%   2  2012-03-15        BryanP      Original (partial) Code
%   1  2012-03-07        BryanP      Pseudo-code

%% ====== Handle Inputs =====
%--- Handle ADP options
if nargin < 2 || isempty(adp_opt)
    adp_opt = struct([]);
end

adp_defaults = {
                %adpSampleBackInd specific settings
                'sbi_state_samples_per_time'            100     % Number of state samples per time period
                'sbi_decisions_per_sample'              10      % Number of decision samples per state
                'sbi_uncertain_samples_per_post'        10     % Number of random/uncertainty samples per time, used for all decisions
                'sbi_min_ok_for_ev'                     0.7    % Min fraction of valid next_predec vals for computing exp_val. setting 0 requires only a single valid val
                %'sbi_build_pre_value_fun'               false  %Explicitly build pre decision value functions for comparision to dpBI. Uses same function approximation as post decision (Not yet implemented)

                %(post decision) value function defaults
                'vfun_approx'           'LocalAvg'
                'vfun_setup_params'     {'AutoExpand', true}
                'vfun_approx_params'    {}

                %decision function defaults
                'decfun_params'         {}
                
                %uncertainty sampling defaults
                'sample_opt'            {'sobol'}
                
                %manage collapsing of non-unique states
                'ops_collapse_duplicate'        true	%If true, reduces state list to unique_states before calling ops. Set to false with stochastic ops costs (e.g. many partially hidden Markov)
                'next_pre_collapse_duplicate'	true	%If true, reduces next-predecision state list to unique_states before computing next period look ahead portion of post decision value function

                %general ADP options used here
                'verbose'               50
                'parallel'              false
                'fix_rand'              false
                'fix_rand_is_done'      false
               };

adp_opt = DefaultFields(adp_opt, adp_defaults);


%--- Handle problem setup
% NOTE: must occur after adp_opts setup, because we use some of those
% settings

% First, announce what were doing
if adp_opt.verbose
    s_num_string = strtrim(sprintf('%g ', adp_opt.sbi_state_samples_per_time/1000));
    fprintf('Sampled Backward Induction ([%s] ksamples/period)\n', s_num_string)
end

% Check required problem fields & fill other defaults
verifyProblemStruct(problem);

%% ====== Standardized Setup =====
% Configure random numbers, including setup consistent stream if desired
adp_opt = utilRandSetup(adp_opt);

% Initialize parallel pool if required (and suppress parallel pool
% initializtion if parallel off)
[cache_par_auto, ps] = utilParSetup(adp_opt);    

% Add additional required functions as needed
problem.fOptimalDecision = utilFunctForProblem(problem, 'fOptimalDecision', @FindOptDecFromVfun);
problem.fRandomSample = utilFunctForProblem(problem, 'fRandomSample', ...
    @(t, n_sample) RandSetSample(problem.random_items, t, n_sample, adp_opt.sample_opt{:}));

%% ====== Additional Setup =====
%Pad samples per period if shorter than number of periods
if length(adp_opt.sbi_state_samples_per_time) < problem.n_periods+1
    adp_opt.sbi_state_samples_per_time(end:(problem.n_periods+1)) = ...
        adp_opt.sbi_state_samples_per_time(end);
end
if length(adp_opt.sbi_decisions_per_sample) < problem.n_periods+1
    adp_opt.sbi_decisions_per_sample(end:(problem.n_periods+1)) = ...
        adp_opt.sbi_decisions_per_sample(end);
end
if length(adp_opt.sbi_uncertain_samples_per_post) < problem.n_periods+1
    adp_opt.sbi_uncertain_samples_per_post(end:(problem.n_periods+1)) = ...
        adp_opt.sbi_uncertain_samples_per_post(end);
end

%-- Initialize POST decision value function (or use the one passed in)
if nargin < 3 || isempty(post_vfun)
    if adp_opt.verbose
        fprintf('    Creating empty post-decision value functions (%s)\n', adp_opt.vfun_approx)
    end

    %Remove "empty post_vfun" to avoid type mismatch issues
    clear('post_vfun')

    %Setup one post decision value function object per time period
    vfun_constructor = str2func(['fa' adp_opt.vfun_approx]);
    for t = 1:problem.n_periods+1
        % Note place holders for initial point and value lists
        post_vfun(t) = vfun_constructor([],[], adp_opt.vfun_setup_params{:});
        
        % Copy names to value fuction 
        post_vfun(t).PtDimNames = problem.state_set{t}.pt_dim_names;
    end
    
    %Flag that we are NOT using old_results (post_vfun) in our options structure
    adp_opt.old_results = false;
else
    if adp_opt.verbose
        fprintf('    Using (copy of) existing post-decision Value Function\n')
    end
    %Make a true, local copy of the value function (since they are handle
    %objects. Otherwise, our additions will effect the stored results)
    for t_idx = problem.n_periods+1:-1:1
        post_vfun(t_idx) = copy(post_vfun(t_idx));
    end

    %Flag that we are using old_results (post_vfun) in our options structure
    adp_opt.old_results = true;
end

%-- Determine the dimensions (length) for states
ndim.pre = zeros(1, problem.n_periods+1);
ndim.post = zeros(1, problem.n_periods+1);
for t = 1:problem.n_periods+1
    ndim.pre(t) = problem.state_set{1}.N_dim;
    %TODO: actually determine post decision size
    ndim.post(t) = ndim.pre(t);
end

%% ====== Sampled Backward Induction Algorithm =====
%% --- Sample terminal states (T+1) ---
t = problem.n_periods + 1;
% TODO: Use problem specific sample function if one is defined
if adp_opt.verbose
    fprintf('    T=%d (terminal period): ', t)
end

% Actually sample states. 
% Note: pre- & post-decision states are the same for the terminal value
pre_state_list = sample(problem.state_set{t}, adp_opt.sbi_state_samples_per_time(t));

% -- Compute contribution for sampled states
% Note: Memoized operations in terminal period must be handled by user. Typically
% by calling a memoized ops function in fTerminalValue and storing the Ops
% table in a handle derived class.
state_values = problem.fTerminalValue(problem.params, t, pre_state_list);

%Add these points to the function approximation
if not(isnan(problem.fMapState2Vfun))
    pre_state_list = problem.fMapState2Vfun(problem, pre_state_list, t);
end
post_vfun(t).update(pre_state_list, state_values);

if adp_opt.verbose
    fprintf('Done\n')
end

%% --- Run earlier periods ---
% Can't easily parallelize this outer loop b/c need future time steps results as we work backward 
%  Note: n_periods is number of decision periods, so final period is n+1.
for t = problem.n_periods:-1:1
    if adp_opt.verbose
        fprintf('    T=%d:', t)
    end

    %Sample pre-decision states
    % will be used for building post_vfun
    pre_state_list = sample(problem.state_set{t}, adp_opt.sbi_state_samples_per_time(t));

    
    %--- Construct post decision value function ---
    
    %-- First obtain a sample of post decision states
    % Since we only have pre-decision states listed, we have to build a set
    % of post-decision states using pre-decision + decisions
    if adp_opt.verbose
        fprintf('S')
    end

    % find decisions for each (in parallel if enabled)
    %
    % First cache relevant portions of large structures to reduce worker
    % communication by not having to send the entire structure 
    fn_decision_set = problem.fDecisionSet;
    fn_decision_apply = problem.fDecisionApply;
    params_only = problem.params;
    dec_per_sample = adp_opt.sbi_decisions_per_sample(t);
    samp_per_time = adp_opt.sbi_state_samples_per_time(t);
    % Pre-allocate output storage. 
    % Use cell since may have multiple entries per state sample and parfor
    % doesn't support the required complex indexing 
    post_state_list = cell(samp_per_time,1);
    decision_list = cell(samp_per_time,1);
    
    parfor samp = 1:samp_per_time
        %Extract valid decisions
        decision_set = fn_decision_set(params_only, t, pre_state_list(samp, :)); %#ok<PFBNS>
        %And sample these choices
        decision_list{samp} = sample(decision_set, dec_per_sample);
        %Build a piece of the post_decision sample
        post_state_list{samp} = fn_decision_apply(params_only, t, pre_state_list(samp, :), decision_list{samp}); %#ok<PFBNS>
    end
    
    % Finally reconstruct the possibly parallel pieces into a full list
    post_state_list = cell2mat(post_state_list);
    decision_list = cell2mat(decision_list);
        
    %-- Compute values for every post_state
    % Explanation: The goal here is to estimate the cost-to-go/value for
    % each sampled post_state. A representative basic algorithm would be:
    %
    %>>>   for each post_state
    %>>>     compute after_decision operations cost
    %>>>   for every sampled post_state (not just unique to maintain sample diversity) 
    %>>>     sample uncertainty and store change
    %>>>     compute uncertainty_contribution
    %>>>     compute after random operations costs
    %>>>     create resulting next_pre_state_list
    %>>>     for each next_pre_state
    %>>>       FindOptimalDecision (based on decision_costs and next post_decision value function) 
    %>>>       compute next_pre_value as sum(optimal_decision_cost + optimal_post_vfun_value) 
    %>>>       compute post_random_val as sum( uncertainty_contribution, after_random_ops, (1-disc_rate) * next_pre_value) 
    %>>>     compute E[random_value] as sum(post_random_val)/n_unique_next_pre 
    %>>>     compute total post_decision_value as sum(after_dec_ops_cost, expected_random_value)
    %>>>   build post decision value function for current time
    %
    % However FindOptimalDecision is typically an expensive computation,
    % and is also fairly likely that there will be overlap in pre_state
    % lists when multiple post_state and uncertainty combinations lead to
    % the same pre_state. As a result the nested loops are split by first
    % computing the unique set of pre_states across all post_states, then
    % finding the optimal decision and future values, and finally matching
    % these future values with the originating post_state to compute the
    % expected values.
    %
    % Also note, that here we don't obtain probabilities from the problem's
    % random set, since we are using Monte Carlo and assume the probability
    % is proportional to the number of samples
    
    % Initialize storage
    % Here we maintain a cell array so that each post_state can track
    % multiple uncertainty/next_decision/next_value combinations
    % Note: also ensures desired cell column vector shape
    n_post_states = size(post_state_list,1);
    next_pre_list = cell(n_post_states,1);
    uncertainty_list = cell(n_post_states,1);
    decision_list_for_pre = cell(n_post_states,1);
    uncertainty_contrib = cell(n_post_states,1);

    % Cache only required structure pieces for use in the following set of
    % parfor loops
    % Note: some values, e.g. params already cached
    fn_random_apply = problem.fRandomApply;
    fn_random_cost = problem.fRandomCost;
    fn_random_sample = problem.fRandomSample;
    fn_optimal_decision = problem.fOptimalDecision;
    %Note: Decision costs not included b/c assume internal loop for computing operations cost (with memoization?) 
    n_periods = problem.n_periods;
    adp_verbose = adp_opt.verbose;
    rand_per_post_state = adp_opt.sbi_uncertain_samples_per_post(t);
    vfun_approx_params = adp_opt.vfun_approx_params;
    this_vfun = post_vfun(t+1);

    %>>>   for each post_state
    % Note: for loop is implicit, since computed internal to fOps*
    
    % Merge unique post states if desired
    if adp_opt.ops_collapse_duplicate
        % Keep only unique values (since multiple copies of same state will
        % be treated the same.
        [post_state_list, unique_post_map, unique_to_full_post_map] = ...
            uniquetol(post_state_list, 'ByRows', true);
        decision_list = decision_list(unique_post_map, :);
    end

    %>>>     compute after_decision operations cost
    if not(isempty(problem.fOpsAfterDecision))
        % Compute unique ops costs, using internal Ops loop
        after_decision_ops = problem.fOpsAfterDecision(params_only, t, ...
            post_state_list, decision_list);
    else
        after_decision_ops = 0;
    end

    % Rebuild full lists if needed
    if adp_opt.ops_collapse_duplicate
        post_state_list = post_state_list(unique_to_full_post_map, :);
        decision_list = decision_list(unique_to_full_post_map, :);
        after_decision_ops = after_decision_ops(unique_to_full_post_map, :);
    end
    
    %>>>   for every sampled post_state (not just unique to maintain sample diversity) 
    %Loop to find full set of next_pre_decision states
    
    % Note: since we typically sample multiple times per post decision
    % state, each result is stored as a cell array
    n_post_state = size(post_state_list, 1);
    parfor post_idx = 1:n_post_state
        if adp_verbose
            DisplayProgress(adp_verbose,post_idx)
        end
        
        %>>>     sample uncertainty and store change
        % Sample Random outcomes to get to next pre-states
        uncertainty_list{post_idx} = fn_random_sample(t, rand_per_post_state); %#ok<PFBNS>, because OK to broadcast reduced fn_* variable
        next_pre_list{post_idx} = fn_random_apply(params_only, t, post_state_list(post_idx,:), uncertainty_list{post_idx});  %#ok<PFBNS>
        decision_list_for_pre{post_idx} = repmat(decision_list(post_idx, :), size(uncertainty_list{post_idx}, 1), 1);

        %>>>     compute uncertainty_contribution
        if not(isempty(fn_random_cost))
            uncertainty_contrib{post_idx} = fn_random_cost(params_only, t, post_state_list(post_idx,:), uncertainty_list{post_idx});
        else
            uncertainty_contrib{post_idx} = 0;
        end
    end
    
    % Concatinate lists of states, decisions, and uncertainty into a single
    % matrix. Can re-split later using num2mat (or some tensor mapping)
    uncertainty_list = cell2mat(uncertainty_list);
    next_pre_list = cell2mat(next_pre_list);
    decision_list_for_pre = cell2mat(decision_list_for_pre);
    
    uncertainty_contrib = cell2mat(uncertainty_contrib);
    
    if adp_opt.next_pre_collapse_duplicate || adp_opt.ops_collapse_duplicate
        % Identify unique next_pre_states, extract corresponding maps and
        % collapse set of next states to explore
        [next_pre_list, pre_map, next_pre_to_post_map] = uniquetol(next_pre_list, 'ByRows', true);
        decision_list_for_pre = decision_list_for_pre(pre_map, :);
        uncertainty_list = uncertainty_list(pre_map, :);
    end
    
    %>>>     compute after_random operations cost
    if not(isempty(problem.fOpsAfterRandom))
        % Compute unique ops costs, using internal Ops loop
        after_random_ops = problem.fOpsAfterRandom(params_only, t, ...
            next_pre_list, decision_list_for_pre, uncertainty_list);
    else
        after_random_ops = 0;
    end
    
    
    % Rebuild full post lists now if collapsed only for ops
    if not(adp_opt.next_pre_collapse_duplicate)
        % Rebuild list of states
        next_pre_list = next_pre_list(next_pre_to_post_map, :);
        
        if adp_opt.ops_collapse_duplicate && not(isscalar(after_random_ops))
            % Rebuild full after random ops costs
            after_random_ops = after_random_ops(next_pre_to_post_map, :);
        end
    end

    
    n_next_pre = size(next_pre_list, 1);    
    next_pre_value = zeros(n_next_pre, 1);    
    
    %>>>     for each next_pre_state
    %Gather contributions for this set of next_pre_states
    parfor next_pre_idx = 1:n_next_pre
        
        %>>>       FindOptimalDecision (based on decision_costs and next post_decision value function) 
        %>>>       compute next_pre_value as sum(optimal_decision_cost + optimal_post_vfun_value) 
        if t == n_periods
            % Note: If looking ahead from the final decision period, there
            % will be no next period decision, only terminal values. 
            next_pre_value(next_pre_idx) = ...
                this_vfun.approx(next_pre_list(next_pre_idx, :), vfun_approx_params{:}); %#ok<PFBNS>      
        else
            next_pre_state = next_pre_list(next_pre_idx, :);
            
            %Find optimal next period decision, associated decision
            %contribution (in t+1 money) and post-decision value
            [~, next_dec_contrib, ~, next_post_val] = ...
                fn_optimal_decision(problem, t+1, next_pre_state, this_vfun, vfun_approx_params); %#ok<PFBNS> 

            %Compute post-decision value function for this time period
            %(store in t+1 money)
            next_pre_value(next_pre_idx) = next_dec_contrib + next_post_val;
                
        end
    end
    
    if adp_opt.verbose && not(0 == mod(adp_opt.sbi_state_samples_per_time(t+1), adp_opt.verbose * 50))
            fprintf('%d\n', adp_opt.sbi_state_samples_per_time(t))
    end
    
    % Rebuild full post lists now if collapsed for next pre
    if adp_opt.next_pre_collapse_duplicate
        % Rebuild next period value list
        next_pre_value = next_pre_value(next_pre_to_post_map, :);
       
        if adp_opt.ops_collapse_duplicate && not(isscalar(after_random_ops))
            % Rebuild full after random ops costs
            after_random_ops = after_random_ops(next_pre_to_post_map, :);
        end
    end
    
    %>>>       compute post_random_val as sum( uncertainty_contribution, after_random_ops, (1-disc_rate) * next_pre_value) 
    after_random_value = uncertainty_contrib + after_random_ops + (1-problem.discount_rate) * next_pre_value;

    % report number of invalid after_random values. These typically
    % represent locations where the value function in the next period is
    % too sparse for accurate estimates using current vfun_approx_params,
    % though could also be invalid uncertainty_contrib or after_random_ops
    if adp_opt.verbose && (nnz(isnan(after_random_value)) > 0)
            fprintf('Ignoring %d/%d invalid entries in expectation of post-decision\n', ...
                nnz(isnan(after_random_value)), length(after_random_value))
    end

    %reshape after_random_value to have one row per post_decision state
    % Note: reshape first fills by column so need to create then transpose
    after_random_value = reshape(after_random_value, [], n_post_state)';
    
    %>>>     compute E[random_value] as sum(post_random_val)/n_unique_next_pre 
    expected_random_value = mean(after_random_value, 2, 'omitnan');
    
    % Flag rows with insufficient valid values to get a good expected value
    %Note using strict > (not >=) allows sbi_min_ok_for_ev to be set to 0 to
    %require only one valid forward sample
    invalid_mask = sum(isnan(after_random_value),2) ./ adp_opt.sbi_uncertain_samples_per_post(t) > adp_opt.sbi_min_ok_for_ev;
    expected_random_value(invalid_mask) = NaN;
    
    %>>>     compute total post_decision_value as sum(after_dec_ops_cost, expected_random_value)
    post_decision_value = after_decision_ops + expected_random_value;

    %>>>   build post decision value function for current time
    % First flag and remove NaNs to ensure a clean function approximation
    valid_post_val_map = not(isnan(post_decision_value));

    
    if any(valid_post_val_map)
        %Add any valid points to the function approximation
        if not(isempty(problem.fMapState2Vfun))
            post_state_list(valid_post_val_map,:) = problem.fMapState2Vfun(problem, post_state_list(valid_post_val_map,:), t);
        end
        post_vfun(t).update(post_state_list(valid_post_val_map,:), post_decision_value(valid_post_val_map,:));
    else
        err_msg = sprintf('No valid new post_vfun points at t=%d. Aborting Sampled Backward Induction', t);
        err_id = 'ADP:SBI:NoValidPts';
        %If no valid new points warn when we are extending an existing
        %value function (which may have values for prior periods)
        if adp_opt.old_results
            warning(err_id, err_msg) %#ok<SPWRN> b/c want to reuse msg
            break
        else
            %or abort if no existing value function since we won't be able
            %to get valid earlier values
            error(err_id, err_msg) %#ok<SPERR> b/c want to reuse msg
        end
    end

end
% Note: Finished building function approximation

%-- Find optimal build for first period (single state)
first_state = problem.state_set{1}.as_array();
if size(first_state, 1) > 1
    warning('adp:MultipleFirstStates', 'Multiple initial states defined, using first in list')
    first_state = first_state(1, :);
end

[results.first_decision, decision_contrib, ~, future_val] = ...
    problem.fOptimalDecision(problem, 1, first_state, post_vfun(1), adp_opt.vfun_approx_params);

results.objective = decision_contrib + (1-problem.discount_rate) * future_val;
results.post_vfun = post_vfun;
results.adp_opt = adp_opt;

if adp_opt.verbose && not(0 == mod(adp_opt.sbi_state_samples_per_time(1), adp_opt.verbose * 50))
    fprintf('%d\n',adp_opt.sbi_state_samples_per_time(1))
end

%% ===== Clean-up =====
%Reset Auto-parallel state
if not(isempty(cache_par_auto))
    % when non-empty, we already created ps
        ps.Pool.AutoCreate = cache_par_auto;
end


end % Main Function
