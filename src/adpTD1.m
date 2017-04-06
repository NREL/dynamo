function results = adpTD1(problem, varargin)

% WARNING: INCOMPLETE NOT UPDATED FOR NEW PROBLEM STRUCTURE

% adpTD1 Temporal Difference Optimization with lamba=1 (double pass)
%
% results = adpTD1(problem, adp_options)
%
% Understanding the four types of states in adpTD1:
%  pre-decision: a description of the system state at the beginning of the
%   time period before a decision has been made. In the ADP toolbox this
%   state only contains the subset of state information required to make a
%   decision.
%  full-state: The ADP toolbox distinguishes a second type of pre-decision
%   state that includes the full set of state dimensions required to
%   simulate the problem, which is a superset of the decision influencing
%   pre-decision state.
%  post-decision: A description of the system state immediately
%   after making a decision, typically the post-decision state has the same
%   form and dimensions as the pre-decision state.
%  vfun: A mapping of the post-decision state that is used for value
%   function approximation. It often contains a sub-set of the
%   post-decision state or a transformed set of alternative dimensions.
%
% State Representation:
%   With the exception of the full-state, all states are represented as
%   numeric row vectors. This allows stacking multiple states to form a 2-D
%   numeric array. Full-states can be any form of valid matlab value which
%   allows maintaining sophisticated Simulink or other simulation data.
%
% Required problem fields:
%   TODO: update for new problem structure
%
% Adp options:
%   See ADP defaults definition. Any of these can be overridden by
%   supplying a string value pair as optional arguements or by passing the
%   corresponding structure or cell array form.


%OLD DOCUMENTATION
%
% Required problem fields (or properties if using problem objects)
%   n_periods
%   discount_rate   inter-decision period discount rate
%
% Common problem fields
%   rand        A struct with subfields of RandProc objects for state and
%                decision to be used during bootstrap iterations
%
% Required problem function handle fields (or methods if object-oriented):
%   fSim            Simulates the impact of random processes on the problem
%                    based on the current postdecision state and/or the
%                    current "full" predecision state. The most recent
%                    decision is included for updating the full state.
%                    Returning the next (simplified) pre-decision state and
%                    updated "full" pre-decision state. Also computes the
%                    non-decision contribution of this transition. During
%                    the final time period, should return the terminal
%                    value. Function signature:
%                      [next_pre, sim_contrib, new_full_state]...
%                        = fSim(problem, t, decision, post_state, full_state);
%
%   fApplyDscn      Applies the optimal decision (see functions below) to
%                    convert from the pre-decision to post-decision states.
%                    Function signature. For generalized use, should be able
%                    to handle a vector of states and return the decision
%                    contribution:
%                      [post_states, dec_contrib] =...
%                           fApplyDscn(problem, pre_state, decision_list, t)
%
%   fFirstState     Return the starting state for the simulation. Function
%                    signature:
%                      [first_pre_state, first_full_state] = fFirstState(problem)
%
% Plus Either:
%   fOptimalDec     A function that finds the optimal decision,
%                    post-decision state and related cost values for a
%                    given pre decision state using the provided value
%                    function. Signature:
%                      [decision, dec_contrib, post_state, forward_val] = ...
%                       fOptimalDec(problem, t, pre_s, vfun, adp_options)
%  or
%   fDecision       A function that lists all possible decisions (as rows
%                    in an arbitrary 2-D array) for the current
%                    pre-decision state. This is used automatically by
%                    FindOptDecFromVfun() to exhaustively search the
%                    decision space. Function signature:
%                      d_list = problem.fDecision(problem, pre_s, t)
%
% In addition the following function fields (or methods) can optionally be
% defined for problem specific behavior
%   fPostToVfun     Takes  a set of post decision states (one state per
%                    row) and converts it into a different set of states
%                    (still one per row) for use with the value function.
%                    This allows using alternative bases (e.g. principle
%                    components) or simply a subset of the dimensions. If
%                    omitted, the all postdecision state dimensions are
%                    used. Function signature:
%                      vfun_state_list = fPostToVfun(problem, post_state_list, t)
%   fFullSim        Simulate all forward paths in their entirety based on
%                    a list of decisions obtained during the forward pass
%                    based on the value function. This is useful in two
%                    cases:
%                       1) With separate sub-model approximation, the
%                       current approximation is used during the forward
%                       pass, the full simulations are run, and then the
%                       sub-model approximation & associated contributions
%                       can be updated simultaneously. This is particularly
%                       helpful for expensive sub-models by enabling better
%                       balancing and non-duplication in parallel sub-model
%                       runs.
%
%                       2) If it is possible to run an approximate,
%                       incremental forward simulation with fSim, but the
%                       full sophisticated model requires simulating the
%                       entire time horizon all at once
%                    The function takes/returns 3d values with one row per
%                    time period, one column per associated dimension (eg
%                    each state dimension), and a 3rd dimension with one
%                    entry per parallel path. Function signature:
%                      [sim_contrib_array, problem] = ...
%                        fFullSim(problem, t, pre_state_3d, decision_3d, sim_contrib_3d, post_state_3d);
%                    For complete simulations over all time, t is specified
%                    as [].
%
%   fDecCompare     Comparision function for decisions used for
%                    convergence checks and backpass abort calculations.
%                    Default: isequaln(). Function signature:
%                      t_or_f = fDecCompare(new_decision, old_decision)
%   fValCompare     Comparision function for value functions used for
%                    convergence checks and backpass abort calculations.
%                    Default: simple internal relative value compare with
%                    divide by zero work around. Function signature:
%                      t_or_f_list = fValCompare(new_val_list, old_val_list, tol)
%                    The associate tolerance is set by the 'tol_bkps_val'
%                    adp option

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  21  2017-03-04 23:43  BryanP      Early rework for new problem structure: abstract common util* functions 
%  20  2016-04-15 02:03  BryanP      Update documentation for required problem parameters
%  19  2012-07-06 16:55  BryanP      Added time to PostToVfun calls
%  18  2012-07-05        BryanP      Integrated "smart" sampled boostrap
%  17  2012-07-03 14:55  BryanP      Renamed adpSBI (was adpSampBackInd)
%  16  2012-06-25 21:35  BryanP      updated for renamed SetDefaultOpts (was SetOptions in adp svn <86)
%  15  2012-06-12 11:05  BryanP      Interface clean up:
%                                     - Switch to varargin based options.
%                                     - Move old_results to options
%                                     - Updated comments
%  14  2012-06-08 12:55  BryanP      New adp defaults: decfun_params
%  13  2012-05-15 01:30  BryanP      Corrected 2nd-to-last contribution in back pass
%  12  2012-05-14 12:40  BryanP      Renamed adpTD1 (from parTD1)
%  11  2012-05-11 17:10  BryanP      Cleaned-up discounting
%  10  2012-05-08 14:05  BryanP      Renamed from PreToPost to ApplyDscn
%   9  2012-05-07 11:22  BryanP      Add dimension names to vfun
%   8  2012-05-05 23:55  BryanP      Add adp to results
%   7  2012-05-05 23:50  BryanP      Added support for consistant random numbers (fix_rand)
%   6  2012-04-22 17:30  BryanP      IT WORKS! Corrected backpass value track/update, add obj & p1_desc to results
%   5  2012-04-21 21:00  BryanP      Reordered indicies with s_path last to enable easier slicing
%   4  2012-04-20 06:45  BryanP      Separate postdec from vfun state spaces
%   3  2012-04-16 17:25  BryanP      First complete version
%   2  2012-04           BryanP      Original (partial) Code
%   1  2012-03           BryanP      Pseudo-code


%===========================%
%%     Initialization       %
%===========================%

%Support passing pre-build structures or cell arrays in for options
if length(varargin) == 1
    varargin = varargin{:};
end

adp_defaults = {
                %adp convergence & time values
                'max_iter'              1000
                'sample_per_iter'       10
                'max_time_sec'          Inf

                %adpTemporalDifference related settings
                'explore_iter'          100     %(DISABLED) Frequency of additional exploration mid-optimization
                'explore_alg'           []

                'seed_vfun_init'        'adpSBI'  %String name of function
                'seed_vfun_params'      {'sbi_samples' 100}
                'vfun_update_in_place'  false   %Flag when using as a sub-algorithm to help initialize or explore 
                
                'tol_bkps_val'          1e-4
                'bkps_use_updated_vfun' false
                'bkps_abort'            true
                'bkps_abort_by_dec'     false
                'converge_bkps'         true

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
              
                'old_results'           []      %Existing results (including value function) to start from
               };

adp_opt = DefaultOpts(varargin, adp_defaults);

%==== Additional Setup
%---- Start timer
parTD_time = tic;

%--- Handle problem setup
% NOTE: must occur after adp_opts setup, because we use some of those
% settings

% First announce what we are doing
if adp_opt.verbose
    localDisplayHeader(adp_opt)
end
    
% Check required problem fields & fill other defaults
verifyProblemStruct(problem);

% Add additional required functions as needed
problem.fOptimalDecision = utilFunctForProblem(problem, 'fOptimalDecision', @FindOptDecFromVfun);
problem.fRandomSample = utilFunctForProblem(problem, 'fRandomSample', ...
    @(t, n_sample) RandSetSample(problem.random_items, t, n_sample, adp_opt.sample_opt{:}));

problem.fDecCompare = utilFunctForProblem(problem, 'fDecCompare', @isequaln);
problem.fValCompare = utilFunctForProblem(problem, 'fValCompare', @utilEqualWithTol);

%% >>>>>>>>>>>>>>> EDIT MARKER <<<<<<<<<<<<<<<<<<


%---- Options and value functions
[adp_opt, vfun, iter_start] = InitVfun(adp_opt, problem);
disc_factor = 1-problem.disc_rate;

%---- Initial resutls and log entries
results.vfun = vfun;
results.adp = adp_opt;
results.is_converged = false;
results.n_iter = iter_start;
results.first_desc = NaN;
results.objective = NaN;

results.log.backpass_abort = zeros(1, problem.n_periods+1);

%==========================%
%%   Seed Value Functions  %
%==========================%
if not(isempty(adp_opt.vfun_init_alg))
    vfun = SeedValueFunctions(vfun, adp_opt, problem, results);
end

%====================%
%%   Main Algorithm  %
%====================%
% Run ADP (TD, lambda=1 /double pass)
if adp_opt.verbose
    fprintf('AdpTD1:')
end

%=== Variable Setup
%---- Setup iteration loop related storage
pre_s_list = NaN(problem.n_periods + 1, problem.dims.pre_state, adp_opt.sample_per_iter);
post_s = NaN(problem.n_periods + 1, problem.dims.post_state, adp_opt.sample_per_iter);
decision = NaN(problem.n_periods + 1, problem.dims.decision, adp_opt.sample_per_iter);

dec_contrib = zeros(adp_opt.sample_per_iter, problem.n_periods + 1);
sim_contrib = zeros(adp_opt.sample_per_iter, problem.n_periods + 1);
fwd_val = zeros(adp_opt.sample_per_iter, problem.n_periods + 1);

%---- Setup initial state & state storage
first_pre_state = problem.state_set{1}.as_array();
if size(first_pre_state, 1)
    warning('AdpState:MultipleFirstStates', 'Multiple States defined for t=1. Using only first state')
    first_pre_state = first_pre_state(1,:);
end

% Ensure convergence flag set, for cases where max_iter = 0
is_converged = false;

%% ===== ADP Outer Iteration Loop =======
for iter = (iter_start+1):(iter_start + adp_opt.max_iter)

    %% ===== FORWARD PASS

    %print out progress indicator dots
    DisplayProgress(ceil(adp_opt.verbose/adp_opt.sample_per_iter), iter, [])

    t_start = 1;
    pre_start = repmat(first_pre_state,adp_opt.sample_per_iter,1);

    start_decision = zeros(adp_opt.sample_per_iter,1);
    start_post = adp_opt.fn.ApplyDscn(problem, pre_start, 0, t_start);

        if adp_opt.parallel
            parfor s_path = 1:adp_opt.sample_per_iter
                [decision(:, :, s_path), dec_contrib(s_path, :), post_s(:, :, s_path),...
                 fwd_val(s_path, :), pre_s_list(:, :, s_path), sim_contrib(s_path, :)] ...
                    = ForwardPassCore(problem, fn, vfun, adp_opt, ...
                        t_start, pre_start(s_path,:), ...
                        in_bootstrap, start_decision(s_path,:), start_post(s_path, :));
            end
        else
            for s_path = 1:adp_opt.sample_per_iter
                [decision(:, :, s_path), dec_contrib(s_path, :), post_s(:, :, s_path),...
                    fwd_val(s_path, :), pre_s_list(:, :, s_path), sim_contrib(s_path, :)] ...
                    = ForwardPassCore(problem, fn, vfun, adp_opt, ...
                        t_start, pre_start(s_path,:), ...
                        in_bootstrap, start_decision(s_path,:), start_post(s_path, :));
            end
        end

        %% ===== BACKWARD PASS
        % Note: terminal period partly treated like any other b/c descion costs
        % set to zero

        %Initialize back-pass abort tracking assuming all paths good.
        backpass_ok = true(adp_opt.sample_per_iter, 1);

        %Loop over all previous periods and update the value functions
        %  Note: can't parallelize at this level, b/c of inter-period
        %  dependancies.
        % TODO: consider larger parallelization over paths (currently
        % vectorized)
        for t = problem.n_periods+1: -1: t_start
            %Convert postdec states for current time period to vfun
            %space (if needed)
            vfun_state = post_s(t, :, :);
            %Reshape required b/c slicing from 3-D
            vfun_state = reshape(vfun_state, [], adp_opt.sample_per_iter)';
            if not(isempty(problem.MapState2Vfun))
                vfun_state = problem.MapState2Vfun(problem, vfun_state, t);
            end
            
            % Reset accumulated value sum if in terminal period
            %
            % TODO: Update handling of Ops functions
            %
            % The terminal and second-to-last periods are special b/c their
            % post-decision state only includes simulation results and not
            % decision contributions.
            % They also share the same simulation results (terminal period
            % values) but match these to different states: the n_periods+1
            % pre-decision state for the terminal period and the n_periods
            % post-decision state for the final decision period.
            %
            % Therefore:
            if t >= problem.n_periods
                % Initialize/reset actual value accumulator for terminal period
                % to avoid double counting sim_contrib
                actual_vals = zeros(adp_opt.sample_per_iter, 1);
                contrib = sim_contrib(:, t);
%                %and then again for 2nd to last 
%             elseif t == problem.n_periods
%                 actual_vals = zeros(adp_opt.sample_per_iter, 1);
%                 contrib = sim_contrib(backpass_ok, t+1);
            else
            % Compute expected contribution for non-terminal periods
            %
            % Now the actual_vals has already been loaded with the
            % next period actual values and we need to compute a
            % contribution that is a function of the next-period decision
            % and next-period simulation. No discounting required b/c the
            % postdecision value functin is stored in t+1 money
                contrib = sim_contrib(backpass_ok, t+1)...
                           + dec_contrib(backpass_ok, t+1);
            end

            % If we opt to use the updated value function approximation on the
            % backpass, first add the observation and then look up these values
            % Note: current period (t) value function stored in t+1
            % dollars, but need to discount actual (future) values since
            % they are in t+2 money
            if adp_opt.bkps_use_updated_vfun

                % Update the results from the forward pass. Bad values will be NaNs 
                val_to_update = contrib + disc_factor * actual_vals;
                %Update OK (or not) mask with new NaNs
                backpass_ok = backpass_ok & not(isnan(val_to_update));

                %Update the value function for good states
                vfun(t).update(vfun_state(backpass_ok, :), ...
                    val_to_update(backpass_ok));
                
                % Now use the updated value function to get revised estimates
                %  for all (good & bad) forward pass states
                actual_vals = vfun(t).approx(vfun_state, adp_opt.vfun_approx_params{:});

            else
                % If not updating the value function, simply use the observed
                % values from the forward pass and then update the value
                % function
                actual_vals(backpass_ok) = ...
                    contrib + disc_factor * actual_vals(backpass_ok);

                %use ok mask to remove any NaNs to preven NaN in update call
                backpass_ok = backpass_ok & not(isnan(actual_vals));
                vfun(t).update(vfun_state(backpass_ok, :), actual_vals(backpass_ok));
            end

            % Handle backpass abort for each simulation path
            if adp_opt.bkps_abort && not(in_bootstrap)
                new_backpass_ok = backpass_ok;
                
                if adp_opt.bkps_abort_by_dec && t <= problem.n_periods
                %If aborting based on changed decisions:
                %Flag any decisions that have changed too much as not OK
                %Note: no decisions in the terminal state: will only abort those on values
                    if adp_opt.parallel
                        parfor s_path = 1:adp_opt.sample_per_iter
                            if backpass_ok(s_path)
                                new_backpass_ok(s_path) = ...
					BackPassDecCheckCore(problem, fn, vfun(t), adp_opt, t, ...
                                        pre_s_list(t, :, s_path), decision(t, :, s_path)); %#ok<PFBNS>
                            end
                        end
                    else
                        for s_path = 1:adp_opt.sample_per_iter
                            if backpass_ok(s_path)
                                new_backpass_ok(s_path) = ...
					BackPassDecCheckCore(problem, fn, vfun(t), adp_opt, t, ...
                                        pre_s_list(t, :, s_path), decision(t, :, s_path));
                            end
                        end
                    end
                else
                    % Otherwise we check our forward pass vs the current
                    % (post update) values
                    if adp_opt.bkps_use_updated_vfun
                        new_vals = actual_vals(backpass_ok);
                    else
                        new_vals = approx(vfun(t), vfun_state(backpass_ok, :), adp_opt.vfun_approx_params{:});
                    end

                    new_backpass_ok(backpass_ok) = adp_opt.fn.ValCompare(new_vals, fwd_val(backpass_ok, t), adp_opt.tol_bkps_val);

                end
                
                results.log.backpass_abort(t) = results.log.backpass_abort(t)...
                    + nnz(backpass_ok) - nnz(new_backpass_ok);
                backpass_ok = new_backpass_ok;

                % If we no longer have any paths to update, fully abort the
                % backpass
                if not(any(backpass_ok))
                    break
                end
            end

        end


    %% ===== END OF ITERATION CALCS & CHECKS =======
    %---- Check for convergence
    % If we haven't converged on decisions for the backpass, assume we
    % haven't converged
    if not(in_bootstrap) && adp_opt.converge_bkps && all(backpass_ok)
        if adp_opt.verbose
            fprintf('BkPs_')
        end
        is_converged = true;
    end

    % Log progress

    if is_converged
        if adp_opt.verbose
            fprintf('Converged (%d)\n', iter)
        end
        break
    else
        %Check for time out
        if toc(parTD_time) > adp_opt.max_time_sec
            if adp_opt.verbose
                fprintf('Timeout\n')
            end
            break
        end
    end
end

if adp_opt.verbose && not(is_converged) && not(isempty(iter)) && iter == adp_opt.max_iter
    fprintf('\nIteration Limit (%d) reached\n', iter)
end


results.vfun = vfun;
results.adp = adp_opt;
results.is_converged = is_converged;
results.n_iter = iter;
[results.first_desc, dec_contrib, ~, cost_to_go] = adp_opt.fn.OptimalDec(problem, 1, first_pre_state, vfun(1), adp_opt);
results.objective = dec_contrib + disc_factor * cost_to_go;

%% ===== Clean-up =====
    %Reset Auto-parallel state
    if not(isempty(cache_par_auto))
        % when non-empty, we already created ps
            ps.Pool.AutoCreate = cache_par_auto;
    end
end

%% ============ Helper Functions =============
%----------------------
%   localDisplayHeader
%----------------------
function localDisplayHeader(adp_opt)
    fprintf('\nRunning Temporal Difference with lambda=1 (aka double pass)\n')
    if adp_opt.parallel
        par_string = 'on';
    else
        par_string = 'off';
    end
    if isempty(adp_opt.old_results) ...
        || not(isfield(adp_opt.old_results, 'n_iter')) ...
        || adp_opt.old_results.n_iter <=0
        restart_string = '';
    else
        restart_string = ' added';
    end

    fprintf('    max %d%s iterations, %d samples/iter\n', ...
                adp_opt.max_iter, restart_string, adp_opt.sample_per_iter)
    fprintf('    max time %d sec, parallel %s\n', ...
                adp_opt.max_time_sec, par_string)
end

%----------------------
%%   InitVfun
%----------------------
function [adp_opt, vfun, iter_start] = InitVfun(adp_opt, problem)
    % If no existing/old results
    if isempty(adp_opt.old_results)
        %--Setup value function from scratch
        vfun_constructor = str2func(['fa' adp_opt.vfun_approx]);
        %Create array of value functions. Reverse loop to pre-allocate size
        for t = problem.n_periods+1:-1:1
            %Note place holders for initial point and value lists
            vfun(t) = vfun_constructor([],[], adp_opt.vfun_setup_params{:});
            %Copy dimension names from state sets
            vfun(t).PtDimNames = problem.state_set{t}.pt_dim_names;
        end

        iter_start = 0;

    else %If existing value funtion exists, use it
        if isfield(adp_opt.old_results, 'n_iter')
            iter_start = adp_opt.old_results.n_iter;
        else
            iter_start = 0;
        end

        if adp_opt.verbose
            if resume
                fprintf('Resuming at iteration #%d\n', iter_start)
            end
            fprintf('Using (copy of) existing Value Function\n')
        end
        %Make a true, local copy of the value function (since they are handle
        %objects. Otherwise, our additions will effect the stored results)
        for t = problem.n_periods+1:-1:1
            vfun(t) = copy(adp_opt.old_results.vfun(t));
        end

        %Simplify old_results to streamline our options structure
        adp_opt.old_results = true;
    end
end

%----------------------
%   Seed Value Functions
%----------------------
function vfun = SeedValueFunctions(vfun, adp_opt, problem, results)
    % If requested Initial Value function approximation using start_algorithm
    start_alg = str2func(adp_opt.seed_vfun_init);

    % Run the starting algorithm
    if adp_opt.verbose
        fprintf('Seeding value functions using %s', adp_opt.seed_vfun_init)
    end
    
    %Configure to use existing value functions
    % Note: structures are passed by value in MATLAB
    results.vfun = vfun;
    adp_opt.old_results = results;
    adp_opt.vfun_update_in_place = true;    %Use existing value functions
    
    % Actually run the startup algorithm
    seed_results = start_alg(problem, adp_opt);
    vfun = seed_results.vfun;

end

%----------------------
%   ForwardPassCore
%----------------------
function [decision, dec_contrib, post_s, fwd_val, pre_s_list, sim_contrib] ...
            = ForwardPassCore(problem, vfun, adp_opt, ...
                              t_start, pre_start, ...
                              in_bootstrap, start_decision, start_post)
% "Insides" of ForwardPass to be shared by for & parfor loops

    %Setup iteration loop related storage
    pre_s_list = NaN(problem.n_periods + 1, problem.dims.pre_state);
    post_s = NaN(problem.n_periods + 1, problem.dims.post_state);
    decision = NaN(problem.n_periods + 1, problem.dims.decision);

    dec_contrib = NaN(1, problem.n_periods + 1);
    sim_contrib = NaN(1, problem.n_periods + 1);
    fwd_val = NaN(1, problem.n_periods + 1);

    full_state = full_start;
    pre_s_list(t_start, :) = pre_start;

    for t = t_start:problem.n_periods
        pre_s = pre_s_list(t, :);

        %Select decision, associated decision contribution (in
        %current t money) and (estimated) post-decision value (in t+1
        %money)
        if in_bootstrap && (t==t_start)
            decision(t, :) = start_decision;
            post_s(t, :) = start_post;
        else
            [decision(t, :), dec_contrib(t), post_s(t, :), fwd_val(t)] ...
                = adp_opt.fn.OptimalDec(problem, t, pre_s, vfun(t), adp_opt);
        end

        % Simulate forward either by
        %  -- realizing stochastic variables to advance from post_state
        %      to the next pre_state
        %   OR
        %  -- using this decision to simulate the full_state (pre)
        %      forward accounting for any stochasticity
        % In both cases, both the pre_state the revised full_state must
        % be returned, although the user is open to supplying
        % placeholder or empty values if they are not required for a
        % specific problem. In addition, the Sim function must compute
        % the simulation portion of the contribution (decision portion
        % of contribution computed with the OptimalDec)
        %
        % Note: the contribution values may be re-computed/updated/
        % over-written if a FullSim function is provided.
        [pre_s_list(t+1, :), sim_contrib(t), full_state]...
            = adp_opt.fn.Sim(problem, t, decision(t, :), post_s(t, :), full_state);

    end

    %% ===== TERMINAL PERIOD (Forward Pass)
    t = problem.n_periods + 1;
    % Map predecision space to post decision space using a zero
    % decision since no decisions to make in terminal period. This also
    % results in on a single possible "post" decision state and
    % corresponding forward pass estimated value so we can assign it
    % directly
    post_s(t, :)...
        = adp_opt.fn.ApplyDscn(problem, pre_s_list(t, :), 0, t);
    %Compute forward values
    if isempty(adp_opt.fn.PostToVfun)
        vfun_state_list = post_s(t, :);
    else
        vfun_state_list = adp_opt.fn.PostToVfun(problem, post_s(t, :), t);
    end

    if not(in_bootstrap)
        fwd_val(t) = vfun(t).approx(vfun_state_list, adp_opt.vfun_approx_params{:});
    end

    % Now simulate to find the actual values
    [~, sim_contrib(t), ~]...
        = adp_opt.fn.Sim(problem, t, 0, post_s(t, :), full_state);
    % Note: no need to set decision cost to zero, b/c it was already
    % initialized that way.

end

%--------------------------
%   BackPassDecCheckCore
%--------------------------
function new_backpass_ok = ...
    BackPassDecCheckCore(problem, vfun, adp_opt, t, pre_s, old_dec)
    new_dec = adp_opt.fn.OptimalDec(problem, t, pre_s, vfun, adp_opt);
    if not(adp_opt.fn.DecCompare(new_dec, old_dec))
        new_backpass_ok(s_path) = false;
    end
end
