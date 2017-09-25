function vector = utilExtendRowVector(vector, n)
% utilExtendRowVector helper function to extend row vectors by replicating the final entry until desired length 
%
% Note: does not truncate if provided vector is already long enough
%
% Examples/doctest
%
% >> utilExtendRowVector([1 2], 3)
% ans =
%      1     2     2
%
% >> utilExtendRowVector([1 2 3], 2)
% ans =
%      1     2     3
%
% >> utilExtendRowVector([1 2 3], 3)
% ans =
%      1     2     3
%
% >> utilExtendRowVector([1 2 3], 4)
% ans =
%      1     2     3     3
% 
% >> utilExtendRowVector([1 2; 10 20], 3)
% ans =
%      1     2     2
%     10    20    20
% 
% >> utilExtendRowVector({2, magic(2)}, 3)
% ans =
%   1×3 cell array
%     [2]    [2×2 double]    [2×2 double]

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  1   2017-09-24 20:36  BryanP      Initial version

cur_len = size(vector, 2);
if cur_len < n
    vector(:, cur_len+1 : n) = repmat( vector(:, cur_len), 1, n - cur_len);
end
