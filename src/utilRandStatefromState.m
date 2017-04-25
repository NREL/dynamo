function cur_rand_state_cell = utilRandStatefromState(full_state, random_state_map_cell)
% UTILRANDSTATEFROMSTATE helper to extract the random portion of a full state for dynamic programming 
% 
% cur_rand_state = utilRandStatefromState(full_state, random_state_map_cell)
%
% Example/tests:
% >> full_state = [1 20 3 44];
% >> utilRandStatefromState(full_state, {[4] [2,1] []})
% 
% ans =
% 
%     [44]    [1×2 double]    [1×0 double]
% 
% >> ans{2}
% 
% ans =
% 
%     20     1
% 
% >> utilRandStatefromState(full_state, {[] [] []})
% 
% ans =
% 
%     [1×0 double]    [1×0 double]    [1×0 double]
%
%
% See also: RandSetSample, RandSetNextJoint, RandProcess

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-04-25 17:05  BryanP      Initial Code

n_rp = length(random_state_map_cell);
cur_rand_state_cell = cell(1, n_rp);
for idx = 1:n_rp
    cur_rand_state_cell{idx} = full_state(1, random_state_map_cell{idx});
end