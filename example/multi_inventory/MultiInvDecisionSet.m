function [ decision_set, params ] = MultiInvDecisionSet(params, t, pre_decision_state) %#ok<INUSL>
%MULTIINVDECISIONSET Produces a set of decisions for the given pre decision state
% Produces a set of decisions for the given pre state (pre_s) at time t
% for use with DP/ADP shared problem structure. 

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   4  2017-04-26 05:35  BryanP      BUGFIX: implement with setCombinWithLimits (finally available)  
%   3  2016-10-27 12:35  BryanP      Added t as second input 
%   2  2016-10-21 15:47  BryanP      Overhaul for new problem structure 
%   1  2010              BryanP      Original version from MIT

% Figure out how much space is left and then use result to build a new set
% for the decisions.
space_remain = params.total_space - (pre_decision_state * params.unit_space');

decision_set = setCombinWithLimits('',  false, params.unit_space, space_remain);

end

