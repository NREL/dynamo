function [state_list, probability_list] = RandSetNextJoint(cell_of_RandProcess, t, cur_rand_state)
% RANDSETJOINT Creates the joint probability from a group of RandProcess objects 
% 
% [state_list, probability_list] = dlistnext (obj, t, state)
%
% See also: RandSetSample, RandProcess

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-04-09 20:05  BryanP      Initial Code 

%cache number of partial distributions to combine
n_rp = length(cell_of_RandProcess);
partial_state_list = cell(1,n_rp);
partial_probability_list = cell(1,n_rp);

%% Divide up current random process states among partials
%extract sizes of random processes to know how to divide cur_rand_state
dim_list = zeros(1, n_rp);
for idx = 1:n_rp
    dim_list(idx) = cell_of_RandProcess{idx}.N_dim;
end

%actually divide up cur_rand_state
cur_rand_state = mat2cell(cur_rand_state, 1, dim_list);

%% Extract partial state and probability lists
for idx = 1:n_rp
    [partial_state_list{idx}, partial_probability_list{idx}] = ...
        cell_of_RandProcess{idx}.dlistnext(t, cur_rand_state{idx});
end

%% Create joint lists by adding one partial at a time to the current joint
% Initialize using the first partial as the joint
state_list = partial_state_list{1};
probability_list = partial_probability_list{1};
% Add remaining joints
for idx = 2:n_rp
    %id sizes for duplicating
    joint_len = size(state_list, 1);
    partial_len = size(partial_state_list{idx}, 1);

    %build state list by duplicating the joint and partial states
    extended_joint_state = repelem(state_list,partial_len, 1);
    extended_partial_state = repmat(partial_state_list{idx}, joint_len, 1);
    state_list = [ extended_joint_state,  extended_partial_state ];

    %build probability list
    % note: save multiplication to end to reduce rounding errors
    extended_joint_probability = repelem(probability_list, partial_len, 1);
    extended_partial_probability = repmat(partial_probability_list{idx}, joint_len, 1);
    probability_list = [extended_joint_probability, extended_partial_probability];
end

% Run multiply to compute resulting probabilities just one time
probability_list = prod(probability_list, 2);
