function [random_outcome_list, probability_list] = RandSetNextJoint(cell_of_RandProcess, t, cur_rand_state)
% RANDSETJOINT Creates the joint probability from a group of RandProcess objects 
% 
% [random_outcome_list, probability_list] = RandSetNextJoint(cell_of_RandProcess, t, cur_rand_state)
%
% See also: RandSetSample, RandProcess

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-04-09 20:05  BryanP      Initial Code 

%% Handle input arguements

% Extract time from first random process if not specified
if nargin < 2 || isempty(t)
    t = cell_of_RandProcess{1}.t;
end
    
% Divide up current random process states among partials
if nargin < 3 
    cur_rand_state = [];
end
%get cur_rand_state as a cell vector
cur_rand_state = utilRandSetCurState(cell_of_RandProcess, cur_rand_state);

%% Extract partial state and probability lists
n_rp = length(cell_of_RandProcess);
partial_state_list = cell(1,n_rp);
partial_probability_list = cell(1,n_rp);

for idx = 1:n_rp
    [partial_state_list{idx}, partial_probability_list{idx}] = ...
        cell_of_RandProcess{idx}.dlistnext(t, cur_rand_state{idx});
end

%% Create joint lists by adding one partial at a time to the current joint
%TODO: compare to using allcomb(). Guessing this might be faster b/c of
%new-ish repelem
%
% Initialize using the first partial as the joint
random_outcome_list = partial_state_list{1};
probability_list = partial_probability_list{1};
% Add remaining joints
for idx = 2:n_rp
    %id sizes for duplicating
    joint_len = size(random_outcome_list, 1);
    partial_len = size(partial_state_list{idx}, 1);

    %build state list by duplicating the joint and partial states
    extended_joint_state = repelem(random_outcome_list,partial_len, 1);
    extended_partial_state = repmat(partial_state_list{idx}, joint_len, 1);
    random_outcome_list = [ extended_joint_state,  extended_partial_state ];

    %build probability list
    % note: save multiplication to end to reduce rounding errors
    extended_joint_probability = repelem(probability_list, partial_len, 1);
    extended_partial_probability = repmat(partial_probability_list{idx}, joint_len, 1);
    probability_list = [extended_joint_probability, extended_partial_probability];
end

% Run multiply to compute resulting probabilities just one time
probability_list = prod(probability_list, 2);
