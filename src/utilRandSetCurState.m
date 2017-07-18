function cur_rand_state = utilRandSetCurState(cell_of_RandProcess, cur_rand_state)
% UTILCRANDSETCURSTATE helper to setup current state cell vector from a RandSet of RandProcess objects 
% 
% cur_rand_state = utilRandSetCurState(cur_rand_state, cell_of_RandProcess)
%
% Example/tests:
% >> rp_set = {rpDiscreteSample({1}), rpDiscreteSample({[10 20; 30 40]})};
% >> utilRandSetCurState(rp_set)
% 
% ans = 
% 
%   1×2 cell array
% 
%     [1]    [1×2 double]
% 
% >> utilRandSetCurState(rp_set, {1, [100, 200]})
% 
% ans = 
% 
%   1×2 cell array
% 
%     [1]    [1×2 double]
% 
% >> utilRandSetCurState(rp_set, [1 100 200])
% 
% ans = 
% 
%   1×2 cell array
% 
%     [1]    [1×2 double]
%
%
% See also: RandSetSample, RandSetNextJoint, RandProcess

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   2  2017-07-14 20:08  BryanP      Update doctests for new output style 
%   1  2017-04-10 17:05  BryanP      Extracted and expanded from RandSetNetJoint v1 

n_rp = length(cell_of_RandProcess);

if nargin < 2 || isempty(cur_rand_state)
    % Extract current states from random process if not given
    cur_rand_state = cell(1, n_rp);
    for idx = 1:n_rp
        cur_rand_state{idx} = cell_of_RandProcess{idx}.cur_state;
    end
else
    % Or if given divide into a cell vector
    if iscell(cur_rand_state)
        %if already a non-empty cell vector assmue its OK
        return
    end
    
    %extract sizes of random processes to know how to divide cur_rand_state
    dim_list = zeros(1, n_rp);
    for idx = 1:n_rp
        dim_list(idx) = cell_of_RandProcess{idx}.N_dim;
    end
    %actually divide up cur_rand_state
    cur_rand_state = mat2cell(cur_rand_state, 1, dim_list);
end
