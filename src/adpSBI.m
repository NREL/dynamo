function results = adpSBI(problem, adp_opt, post_vfun, fn)
% adpSBI Sampled Backward Induction for estimated ADP and intialization
%
% results = adpSBI(problem, adp_opt, post_vfun, fn)
%
% Note: This method is currently only applicable to problem classes where it is
% possible to sample the post decision state without requiring a full,
% simulation derived state (ie it does not work for partially hidden
% Markov/quasi-Markov processes).
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
                'sbi_state_samples_per_time'            20     % Number of state samples per time period
                'sbi_decisions_per_sample'              3      % Number of decision samples per state
                'sbi_uncertain_samples_per_post'        10      % Number of random/uncertainty samples per time, used for all decisions
                'sbi_min_ok_for_ev'                     0.7    % Min fraction of valid next_predec vals for computing exp_val. setting 0 requires only a single valid val

                %value function defaults
                'vfun_approx'           'LocalRegr'
                'vfun_setup_params'     {'AutoExpand', true}
                'vfun_approx_params'    {}

                'decfun_params'         {}
                
                %uncertainty sampling defaults
                'sample_opt'            {}
                
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

% Check required problem fields & fill other defaults
verifyProblemStruct(problem);

% Add additional required functions as needed
if isempty(problem.fOptimalDecision)
    problem.fOptimalDecision = @FindOptDecFromVfun;
end
if isempty(problem.fRandomSample)
    problem.fRandomSample = @(t, n_sample) RandSetSample(problem.random_items, t, n_sample, adp_opt.sample_opt{:});
end

%% ====== Additional Setup =====
% Announce what were doing
if adp_opt.verbose
    s_num_string = strtrim(sprintf('%g ', adp_opt.sbi_state_samples_per_time/1000));
    fprintf('Sampled Backward Induction ([%s] ksamples/period)\n', s_num_string)
end

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

% Setup consistent random numbers if desired
if adp_opt.fix_rand && not(adp_opt.fix_rand_is_done)
    disp('  RAND RESET: adpSBI using fixed random number stream')
    % Reset the random number generator
    rng('default');

    %Old (pre-R2011a) syntax
	%rand_num_sequence = RandStream('mt19937ar');
	%RandStream.setDefaultStream(rand_num_sequence)

    %Indicate that we have already reset the random number generator to
    %prevent sub-functions from resetting again
    adp_opt.fix_rand_is_done = true;
end

% Initialize parallel pool if required (and suppress parallel pool
% initializtion if parallel off)
%
% TODO: figure out which version of MATLAB is required for AutoCreate and
% add alternate code if needed
% Check for parallel programming toolbox (works for R2016a)
cache_par_auto = [];
if exist('parallel', 'dir')
    try
        ps = parallel.Settings;
        % In current (tested with R2016a) versions of MATLAB, calling
        % parfor will automatically initialize a local worker pool, unless
        % the AutoCreate option is disabled
        cache_par_auto = ps.Pool.AutoCreate;
        ps.Pool.AutoCreate = false;
        p = gcp('nocreate');
        if adp_opt.parallel
            if isempty(p)
                parpool();
            else
                if adp_opt.verbose
                    fprintf(' Using existing Parallel Pool\n')
                end
            end
        else
            if not(isempty(p))
                warning('Adp:NotDisablingParallel', 'A Parallel Pool exists so adpSBI will run in parallel despite adp option for non-parallel run')
            end
        end

    catch
        warning('adpSBI:ParallelInitError', 'Error initializing parallel settings')
        lasterr %#ok<LERR>
    end
else
    if adp_opt.parallel
        warning('adpSBI:NoParallel', 'Cannot find parallel setup. Parallel execution will not be used')
    end
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
        %Finally build a piece of the post_decision sample
        post_state_list{samp} = fn_decision_apply(params_only, t, pre_state_list(samp, :), decision_list{samp}); %#ok<PFBNS>
    end
    
    % Finally reconstruct the possibly parallel pieces into a full list
    post_state_list = cell2mat(post_state_list);
    decision_list = cell2mat(decision_list);
    % And keep only unique values (since multiple copies of same state will
    % be treated the same.
    [post_state_list, post_map] = uniquetol(post_state_list, 'ByRows', true);
    decision_list = decision_list(post_map, :);
        
    %-- Compute values for every post_state
    % Explanation: The goal here is to estimate the cost-to-go/value for
    % each sampled post_state. A representative basic algorithm would be:
    %
    %   for each post_state
    %     sample uncertainty and store change AND probability
    %     compute uncertainty_contribution
    %     compute after random operations costs
    %     compute resulting next_pre_state_list
    %     for each next_pre_state
    %       FindOptimalDecision (based on decision_costs and post_decision value function) 
    %       compute value as sum(uncertain_contribution + optimal_decision_cost + optimal_post_vfun_value 
    %     compute expected value
    %
    % However FindOptimalDecision is typically an expensive computation,
    % and is also fairly likely that there will be overlap in pre_state
    % lists when multiple post_state and uncertainty combinations lead to
    % the same pre_state. As a result the nested loops are split by first
    % computing the unique set of pre_states across all post_states, then
    % finding the optimal decision and future values, and finally matching
    % these future values with the originating post_state to compute the
    % expected values.
    
    % Initialize storage
    % Here we maintain a cell array so that each post_state can track
    % multiple uncertainty/next_decision/next_value combinations
    n_post_states = size(post_state_list,1);
    next_pre_list = cell(n_post_states,1);
    uncertainty_list = cell(n_post_states,1);
    decision_list_for_pre = cell(n_post_states,1);
    val_list = cell(n_post_states,1);
    exp_vals = zeros(n_post_states,1);

    % Cache only required structure pieces for use in the following set of
    % parfor loops
    % Note: some values, e.g. params already cached
    fn_random_apply = problem.fRandomApply;
    fn_random_cost = problem.fRandomCost;
    fn_random_sample = problem.fRandomSample;
    fn_ops_after_random = problem.fOpsAfterRandom;
    fn_ops_before_decision = problem.fOpsBeforeDecision;
    n_periods = problem.n_periods;
    adp_verbose = adp_opt.verbose;
    rand_per_post_state = adp_opt.sbi_uncertain_samples_per_post(t);
    vfun_approx_params = adp_opt.vfun_approx_params;
    this_vfun = post_vfun(t+1);

    %Loop to find full set of next_pre_decision states
    parfor post_idx = 1:size(post_state_list, 1)
        if adp_verbose
            DisplayProgress(adp_verbose,post_idx)
        end
        
        % Sample Random outcomes to get to next pre-states
        uncertainty_list{post_idx} = fn_random_sample(t, rand_per_post_state); %#ok<PFBNS>
        next_pre_list{post_idx} = fn_random_apply(params_only, t, post_state_list(post_idx,:), uncertainty_list{post_idx});  %#ok<PFBNS>
        val_list{post_idx} = fn_random_cost(params_only, t, post_state_list(post_idx,:), uncertainty_list{post_idx});  %#ok<PFBNS>
        decision_list_for_pre{post_idx} = repmat(decision_list(post_idx, :), size(uncertainty_list{post_idx}, 1), 1);
    end

    % Identify unique next_pre_states
    next_pre_list = cell2mat(next_pre_list);
    [next_pre_list, pre_map, next_pre_to_post_map] = uniquetol(next_pre_list, 'ByRows', true);
    decision_list_for_pre = cell2mat(decision_list_for_pre);
    decision_list_for_pre = decision_list_for_pre(pre_map, :);
    uncertainty_list_for_pre = cell2mat(uncertainty_list);
    uncertainty_list_for_pre = uncertainty_list_for_pre(pre_map, :);
    
    n_next_pre = size(next_pre_list, 1);    
    val_next_pre_list = zeros(n_next_pre, 1);    
    
    %Compute all types of contributions for this set of next_pre_states
    for pre_idx = 1:n_next_pre
        %Compute operation cost for existing time period after random 
        if not(isempty(fn_ops_after_random))
            val_next_pre_list(pre_idx) = fn_ops_after_random(params_only, t, ...
                next_pre_list(pre_idx, :), decision_list_for_pre(pre_idx, :), uncertainty_list_for_pre(pre_idx, :));
        end
        
        if t==n_periods %Final decision period
            % When looking ahead from this period, there will be no next
            % period decision, only terminal values. So we can directly use
            % the next_pre_list to find values
            %Find values from the next period value function
            val_next_pre_list(pre_idx) = val_next_pre_list(pre_idx) + ...
                this_vfun.approx(next_pre_list(pre_idx, :), vfun_approx_params{:}); %#ok<PFBNS>
        else
%% >>>>>>>>>>>>>>> EDIT MARKER <<<<<<<<<<<<<<<<<<
    %TODO:
    %   Determine (apparent) optimal decision
    %   Compute decision cost and E(future value)
%             %Periods before the final decision period (includes neither the
%             %terminal period (n_periods+1) nor the final decision period
%             %(n_periods)
%             val_next_pre_list(pre_idx) = NaN(size(probability_list));
%             for next_pre_idx = 1:n_next
%                 next_pre_list = next_pre_list(next_pre_idx, :);
%                 %Find optimal next period decision , associated decision
%                 %contribution (in t+1 money) and post-decision value (in
%                 %t+2 money)
%                 [this_desc, next_dec_contrib, this_next, next_post_val] = ...
%                     fn.OptimalDec(problem, t+1, next_pre_list, post_vfun, adp);
% 
%                 %TODO: clever allocation?
%                 if not(isempty(this_desc))
%                     %opt_desc(next_pre_idx,:) = this_desc;
%                     %next_post(next_pre_idx,:) = this_next;
%                     if isempty(fn.FullSim)
%                         %And corresponding non-decision contribution (in T+1 money)
%                         [~, next_sim_contrib] = fn.Sim(problem, t+1, this_desc, this_next, next_pre_list);
%                     else
%                         next_sim_contrib = 0;
%                     end
%                     %Compute post-decision value function for this time period
%                     %(store in t+1 money)
% %                     val_list{state_idx}(next_pre_idx) = next_sim_contrib ...
% %                                           + next_dec_contrib ...
% %                                           + (1-problem.disc_rate) * next_post_val;
%                 end
%             end
        end
        
    % For each pre-state
    %    Select optimial next decision (t+1) 
    %    Compute decision cost and
    %     cost to go based on next post_vfun

%             parfor s = 1:adp_opt.sbi_state_samples_per_time(t);
%                 [vals{s}, prob{s}] = ...
%                     BackIndCore(problem, fn, this_vf, s, t, post_state_list(s,:), adp_opt);
%             end
       
    end
    


    if adp_opt.verbose && not(0 == mod(adp_opt.sbi_state_samples_per_time(t+1), adp_opt.verbose * 50))
            fprintf('%d\n', adp_opt.sbi_state_samples_per_time(t))
    end


        next_opt_dec = vertcat(next_opt_dec{:});
        next_opt_dec = next_opt_dec(pre_map, :);

        next_post = vertcat(next_post{:});
        next_post = next_post(pre_map, :);

        if adp_opt.verbose
            fprintf('Done\n')
        end

        %Simulate contributions for these states
        [sim_contribs, problem] = CPDpFullSim(problem, t+1, next_pre_list, next_opt_dec, [], next_post, adp_opt);

        %Extract the appropriate values for each post state
        for s = 1:adp_opt.sbi_state_samples_per_time(t);
            [~, s_map] = intersect(next_pre_list, next_pre_by_samp{s}, 'rows');
            val_list{s} = val_list{s} + sim_contribs(s_map, :);
        end

  

    %Compute EVs
    %Find possible next pre-descision states for each post state
    for s = 1:adp_opt.sbi_state_samples_per_time(t);
        %Flag NaNs for removal
        valid_val = not(isnan(val_list{s}));

        % Check sample quality
        %Note using strict > allows sbi_min_ok_for_ev to be set to 0 to
        %require only one valid forward sample
        if nnz(valid_val)/size(val_list{s}, 1) > adp_opt.sbi_min_ok_for_ev
            exp_vals(s) = val_list{s}(valid_val)' * probability_list{s}(valid_val);
        else
            exp_vals(s) = NaN;
        end
    end

    %Remove NaNs
    valid_ev = not(isnan(exp_vals));

    % Check sample quality
    if any(valid_ev)
        %Add any valid points to the function approximation
        if not(isempty(fn.PostToVfun))
            post_state_list = fn.PostToVfun(problem, post_state_list(valid_ev,:), t);
        end
        post_vfun(t).update(post_state_list, exp_vals(valid_ev,:));
    else
        err_msg = sprintf('No valid new post_vfun points at t=%d. Aborting Sampled Backward Induction', t);
        err_id = 'ADP:SBI:NoValidPts';
        %Warn & abort if no valid new points
        if nargout > 1
            error(err_id, err_msg) %#ok<SPERR> b/c want to reuse msg
        else
            warning(err_id, err_msg) %#ok<SPWRN> b/c want to reuse msg
            break
        end
    end

end


%-- Finished building function approximation

% If full results requested, find the optimial first period build and build
% the results structure
if nargout > 1
    %-- Find optimal build for first period (single state)
    s = fn.FirstState(problem);

    [results.first_desc, desn_contrib, ~, future_val] = fn.OptimalDec(problem, 1, s, post_vfun(1), adp_opt);

    results.objective = desn_contrib + (1-problem.disc_rate) * future_val;
    results.post_vfun = post_vfun;
    results.adp = adp_opt;
end

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


%% ============ Helper Functions =============
%------------------
%   BackIndCore
%------------------
function [vals, prob, next_pre_list, opt_desc, next_post] = ...
            BackIndCore(problem, fn, post_vfun, s_idx, t, this_post, adp)
% Note: This function allows us to share the "insides" of the forward
% pass loop for both parallel (parfor) and non-parallel (for) uses

end
