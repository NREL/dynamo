function [v, params] = MultiInvTerminalValue(params, t, s_list) %#ok<INUSL>
% MULTIINVTERMINALVALUE DP terminal value for multi-product inventory problem
%
% as described in DP for fFinalValue
%           returns the final value for state s
%
% see also:
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   7  2016-11-10 11:45  BryanP      Renamed write out "Terminal" 
%   6  2016-10-27 12:35  BryanP      Reordered parameters to standardize on t as second input 
%   5  2016-10-21 15:40  BryanP      Reordered input parameters to put params first (for possible future objectifying) 
%   4  2016-04-29 22:11  BryanP      Convert to row vector for each state and support vectors of states 
%   3  2010-08-08 00:30  BryanP      adapted for DP v12 with param returns
%   2  2010-05-15 19:45  BryanP      updated to use params.states
%   1  2010-05-15 09:47  BryanP      adapted from InvTermValue v2

% ------ Compute costs ------
v = s_list * params.term_unit_val';
