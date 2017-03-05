function [cache_par_auto, ps] = utilParSetup(adp_opt)
% UTILPARSETUP Setup parallel workers (or disable) as desired
%
%   [cache_par_auto, ps] = utilParSetup(adp_opt)
%
% Initialize parallel pool if required (and suppress parallel pool
% initializtion if parallel off)
%
% Returns parameters needed to clean up code and return parallel state 
% using something like:
%     %Reset Auto-parallel state
%     if not(isempty(cache_par_auto))
%         % when non-empty, we already created ps
%             ps.Pool.AutoCreate = cache_par_auto;
%     end
% 
% For use with the adp* family of functions.
%
% TODO: figure out which version of MATLAB is required for AutoCreate and
% add alternate code if needed
% Check for parallel programming toolbox (works for R2016a)

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-03-04 23:13  BryanP      Adapted from adpSBI1 v23


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