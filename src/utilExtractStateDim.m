function ndim = utilExtractStateDim(problem)
% UTILEXTRACTSTATEDIM helper function to extract the pre and post state dimensions for adp algorithms

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  1   2017-06-15 04:48  BryanP      Extracted from adpSBI v31 


ndim.pre = zeros(1, problem.n_periods+1);
ndim.post = zeros(1, problem.n_periods+1);
for t = 1:problem.n_periods+1
    ndim.pre(t) = problem.state_set{1}.N_dim;
    %TODO: actually determine post decision size
    ndim.post(t) = ndim.pre(t);
end
