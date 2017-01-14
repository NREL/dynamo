function [s, params] = MultiInvState(i, t, params)
% MULTIINVSTATE DP state iterator function for multi-product inventory problem (modular DP)
%
% as described in DP for fStateIterator...
%           returns a state variable (of any type), for the positive 
%           decimal state index "i"
%                i = 0 returns the number of possible states
%                otherwise the range of i is 1:number of states
%
% The state itself is represented as a row vector of inventory levels for
% each product
%
% see also:
%   DP, MultiInvDecision, MultiInvTermValue, MultiInvCost, MultiInvTransProb, MultiInvInit
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-15 09:47  BryanP      adapted from InvState v2
%   2  2010-05-15 19:45  BryanP      updated to use params.states
%   3  2010-05-22 09:50  BryanP      added support for time-aware states
%   4  2010-08-08 00:30  BryanP      adapted for DP v12 with param returns
%   5  2016-04-28 00:40  BryanP      Added 'all' for ADP use
%   6  2016-04-29 21:40  BryanP      Convert to row vector for each state and support vectors of i 

% Handle special cases:
%  0 to return the number of states
%  'all' to list all states
if isnumeric(i) && length(i) == 1 && i==0
    s = size(params.states,1);
    return
end
if not(isnumeric(i))
    if strcmpi(i, 'all')
        s = params.states;
        return
    else
        error('ADP:MultiInv:InvalidIdx','MultiInv: index must be numeric or ''all''')
    end
end

% return the state corresponding to the index. For this problem the state
% does not chang size with time
s = params.states(i,:);