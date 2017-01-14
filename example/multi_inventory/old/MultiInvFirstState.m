function [first_pre_state, first_full_state] = MultiInvFirstState(problem)
% Simply returns the initial simulation state.  Complies with the "fFirstState"
% function signature required by adpTD1

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2016-04-15 01:03  BryanP      Initial Code
%   2  2016-05-01 04:33  BryanP      Adpted for vector (rather than index) states 

first_pre_state = problem.first_state;
first_full_state = first_pre_state;