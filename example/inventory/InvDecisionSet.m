function [ decision_set, params ] = InvDecisionSet(params, t, current_stock) %#ok<INUSL>
%INVDECISIONSET Produces a set of decisions for the given pre decision state
% Produces a set of decisions for the given pre state (pre_s) at time t
% for use with DP/ADP shared problem structure. 

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  1   2017-09-26 05:35  BryanP      BUGFIX: implement with setCombinWithLimits (finally available)  
%   4  2017-04-26 05:35  BryanP      BUGFIX: implement with setCombinWithLimits (finally available)  
%   3  2016-10-27 12:35  BryanP      Added t as second input 
%   2  2016-10-21 15:47  BryanP      Overhaul for new problem structure 
%   1  2010              BryanP      Original version from MIT

% Figure out how much space is left and then use result to build a new set
% for the decisions.
max_order = params.max_inv - current_stock;

decision_set = setBasic([0; max_order], [], 1);

end

