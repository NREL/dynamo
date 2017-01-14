function d_list = MultiInvDecisionList(params, pre_s, ~)
% Produces a list of decisions for the given pre state (pre_s) at time t
% for use with ADP. Complies with the "fDecision"
% function signature required by adpTD1
%
% Note: Compared to the DP() fTransProb signature, this does not
% (currently) support passing in future_values, and does not easily allow
% updating the problem (as an output)

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2016-04-14 11:50  BryanP      Extracted from MultiInvDecision v6
%   2  2016-04-29 22:11  BryanP      Convert to row vector for each state and support vectors of states 
%   3  2016-05-01 02:31  BryanP      Bug fix: remaining space calculation 

% Figure out how much space is left and then populate state list
space_remain = params.total_space - (pre_s * params.prod_space');
d_list = CombinWithLimits(params.prod_space, space_remain);

