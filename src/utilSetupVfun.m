function [post_vfun, adp_opt] = utilSetupVfun(problem, adp_opt, post_vfun)
% UTILSETUPVFUN helper function to initalize value functions for adp algorithms
%  also displays standarized user notifications

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  1   2017-06-14 07:13  BryanP      Extracted from adpSBI v30 

if nargin<3 || isempty(post_vfun)
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