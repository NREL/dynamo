function [ops_contribs, new_full_states] ...
    = MultiInvOps(params, t, state_list, decision_list, uncertainty_list, full_state_list) %#ok<INUSD,INUSL> extra params listed by name for clarity
% Operations (holding) costs for the multi inventory problem
%
% Implements the shared problem structure function signature for operations
% that could be called in any one (or more) places in the DP/ADP algorithm:
%
%   Type/Timing           State         Decision    Uncertainty
%   ---------------     -------------   --------    -----------
%   Before Decision     pre-decision	empty       empty
%   After Decision      post-decision   provided    empty
%   After Random        next-pre        provided    provided
%
%   Note: must handle (column) vectors for each parameter and each
%   parameter is typically a (row) vector of values

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   3  2016-10-27 13:36  BryanP      Bugfixes for vector inputs, etc. based on doctests 
%   2  2016-10-27 11:36  BryanP      Renamed and overhauled for shared problem structure     
%   1  2016-04-14 11:23  BryanP      Initial Code that uses existing MultiInvTransProb v4 

if nargin < 4 || isempty(decision_list)
    % Before Decision -- nothing to do
    warning('DpMultiInv:UnusedOps','Unused Operations call before decision. Only after decision operations used in Multi-Inventory')
    ops_contribs = zeros(size(state_list,1), 1);
elseif nargin < 5 || isempty(uncertainty_list)
    % After Decision -- now we can compute holding costs
    ops_contribs = state_list * params.hold_cost';
else
    %After Random -- nothing to do
    warning('DpMultiInv:UnusedOps','Unused Operations after random. Only after decision operations used in Multi-Inventory')
    ops_contribs = zeros(size(state_list,1), 1);
end

if nargout > 1
    % Here we would run a more detailed simulation if needed, such as for a
    % partially hidden Markov process
    new_full_states = state_list;
end
