function u_list = RandSetSample(cell_of_RandProcess, n_samples, t, opt)
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
%   3  2017-04-09 20:05  BryanP      Match revised sample() order in our signature, too 
%   2  2017-04-06 00:05  BryanP      Adjusted order of t & n for revised sample
%   1  2016-10-27 14:05  BryanP      Initial Code

if nargin < 3 || isempty(n_samples)
    n_samples = 1;
end

if nargin < 4
    opt = [];
end


% Note: not preallocating b/c not sure of type
u_list = [];
for idx = 1:length(cell_of_RandProcess)
    u_list = [u_list, cell_of_RandProcess{idx}.sample(n_samples, t, opt)]; %#ok<AGROW>
end