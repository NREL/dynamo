function [v, params] = SimpleStoragePvTerminalValue(params, t, s_list) %#ok<INUSL>
% SIMPLESTORAGEPVTERMINALVALUE DP terminal value for simple PV plus storage problem
%
% as described in DP for fFinalValue
%           returns the final value for state s
%
% Notes:
%  -- Calls operations cost for terminal period pre-decision state
%
% see also:
%
% originally by Bryan Palmintier 2017

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-07-16 21:57  BryanP      adapted from MultiInvTerminalValue v7

% ------ Compute costs ------
v = SimpleStoragePvOps(params, t, s_list);
