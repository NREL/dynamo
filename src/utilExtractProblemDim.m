function ndim = utilExtractProblemDim(problem, t)
% UTILEXTRACTPROBLEMDIM helper function to extract problem dimensions (pre, post, decision) for adp algorithms

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  2   2017-06-17 05:03  BryanP      Overhaul to add post and decision sizing and convert all to a single time 
%  1   2017-06-15 04:48  BryanP      Extracted from adpSBI v31 

if nargin < 2 || isempty(t)
    t=1;
end

ndim.pre = problem.state_set{t}.N_dim;
example_pre = problem.state_set{t}.sample();
example_decision_list = problem.fDecisionSet(problem.params, t, example_pre);
example_decision = example_decision_list.sample();
ndim.decision = length(example_decision);
example_post = problem.fDecisionApply(problem.params, t, example_pre, example_decision);
ndim.post = length(example_post);
