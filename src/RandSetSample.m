function u_list = RandSetSample(cell_of_RandProcess, n_samples, t, cur_rand_state)
% RANDSETSAMPLE Draws the specified number of samples from a group of RandProcess objects 
% 
% Produces a row vector for each sample set
%
% Note: no probability returned, since this is a pure monte carlo sample
%
% See also: RandSetNextJoint, RandProcess

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   5  2017-07-14 21:45  BryanP      Remove sample option support 
%   4  2017-04-10 17:55  BryanP      Add support for state dependant samples 
%   3  2017-04-09 20:05  BryanP      Match revised sample() order in our signature, too 
%   2  2017-04-06 00:05  BryanP      Adjusted order of t & n for revised sample
%   1  2016-10-27 14:05  BryanP      Initial Code

%% Handle input arguements

% Assume 1 sample if not specified
if nargin < 2 || isempty(n_samples)
    n_samples = 1;
end

% Extract time from first random process if not specified
if nargin < 3 || isempty(t)
    t = cell_of_RandProcess{1}.t;
end
    
% Divide up current random process states among partials
if nargin < 4 
    cur_rand_state = [];
end
%get cur_rand_state as a cell vector
% JPB changed parameter order
cur_rand_state = utilRandSetCurState(cell_of_RandProcess, cur_rand_state);

%%
% Note: not preallocating b/c not sure of type
u_list = [];
for idx = 1:length(cell_of_RandProcess)
    u_list = [u_list, cell_of_RandProcess{idx}.sample(n_samples, t, cur_rand_state{idx})]; %#ok<AGROW>
end