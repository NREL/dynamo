function [post_state_list] = ...
    MultiInvDecisionApply(params, t, pre_state_list, decision_list) %#ok<INUSL>
% MULTIINVDECISIONAPPLY Applies the provided decision(s) to the current (pre-decision) state 
% 
% 
% Produces a (set of) new post decision state(s).  Complies with
% the "fApplyDscn" function signature required by adpTD1

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   7  2016-10-27 12:31  BryanP      Reordered parameters to standardize on t as second input 
%   6  2016-10-21 15:42  BryanP      Minor tweaks (comments only) for new shared problem structure deployment 
%   5  2016-07-07 00:30  BryanP      Renamed MultiInvApplyDecision to MultiInvDecisionApply for shared problem structure 
%   4  2016-05-01 04:20  BryanP      Allow one state multi-decision, and one decision, mulit-state through bsxfun 
%   3  2016-05-01 03:50  BryanP      Added support for "zero" state for terminal values in adpTD1 
%   2  2016-04-29 21:40  BryanP      Convert to row vector for each state and support vectors of states 
%   1  2016-04-14 11:23  BryanP      Initial Code

% Handle "zero" decision case for adpTD1
if isequal(decision_list, 0)
    post_state_list = pre_state_list;
    return
end

% Find the corresponding next states
post_state_list = bsxfun(@plus, pre_state_list, decision_list);

