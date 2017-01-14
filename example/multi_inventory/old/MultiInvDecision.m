function [d, params] = MultiInvDecision(s, i, t, future_values, params) %#ok<INUSL>
% INVDECISION DP decision iterator function for multi-product inventory problem 
%
% as described in DP for fDecisionIterator...
%           returns a decision, represented by a variable of any type, that 
%           represents one of the possible decisions to be made at state s. 
%                i = '0 returns the number of possible decisions for the 
%                current state. Otherwize the range of i is 1:num decisions
%                for this state
%
% see also:
%   DP, MultiInvState, MultiInvTermValue, MultiInvCost, MultiInvTransProb, 
%   MultiInvInit, CombinWithLimits
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-15 09:47  BryanP      adapted from MultiInvDecision v2
%   2  2010-05-15 19:45  BryanP      converted to use of persistent list
%   3  2010-05-26 07:25  BryanP      updated for ValidCombin v4
%   4  2010-06-07 00:30  BryanP      Allow param update with each sub-function call 
%   5  2010-08-07 11:50  BryanP      adapted for DP v12 with future value functions
%   6  2011-03-29 10:10  BryanP      Renamed MultiInvValidCombin to CombinWithLimits
%   7  2016-04-14 12:14  BryanP      Extract decision list creation to MultiInvDecisionList 
%   8  2016-04-29 22:10  BryanP      Convert to row vector for each state and support vectors of states 

%By using persistent variables, we can keep the values of variables around
%between calls to the function. Here we use this to avoid having to
%recompute the set of valid decisions each time.
persistent decision_list cur_s

% Here we take advantage of the fact that DP algorithms will typically 
% iterate through all of the decisions for a given state sequentially and 
% therefore, we store the list of valid decisions in our persistant
% variables whenever the state changes.
if isempty(cur_s) || not(isequal(s, cur_s))
    cur_s = s;
    decision_list = MultiInvDecisionList(params, s, t);
end    

if i == 0
    d = size(decision_list,1);
    return
end

d = decision_list(i,:);    