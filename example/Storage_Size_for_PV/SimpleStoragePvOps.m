function [ops_values, new_full_states] ...
    = SimpleStoragePvOps(params, t, state_list, decision_list, uncertainty_list, full_state_list) %#ok<INUSD> extra params listed by name for clarity
% SIMPLESTORAGEPVOPS Operations costs for simple storage + PV problem
%
% For the PV + Storage problem, the operations themselves only depend on the
% current state, but are split with user selected scaling for pre & post
% decision states. The terminal (pre-decision) state is scaled
% independently
%
% In all cases the actual operations technical parameters are computed by
% (memoized) calls to OpenDSS and then translated to values here.
%
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
%   2  2017-07-17 17:34  BryanP      Filled in cost computation portion. Works with lookup table 
%   1  2017-07-17 08:43  BryanP      Adapted MultiInvOps v3 

%-1- Start by selecting the appropriate scaling:
if nargin < 4 || isempty(decision_list)
    % Before Decision -- simulate with pre-decision state
    if t == params.n_periods + 1 %Terminal period
        ops_scale = params.ops_terminal_yrs;
    elseif t == 1
        % No pre-decision ops for initial period
        ops_values = zeros(size(state_list,1), 1);
        new_full_states = state_list;
        return
    else
        ops_scale = params.ops_pre_decision_yrs;
    end
elseif nargin < 5 || isempty(uncertainty_list)
    % After Decision -- simulate with post-decision state
    ops_scale = params.ops_post_decision_yrs;
else
    %After Random -- nothing to do
    warning('DpSimpleStoragePv:UnusedOps','Operations after random not used in Simple Storage + PV')
    ops_values = zeros(size(state_list,1), 1);
    new_full_states = state_list;
    return
end

%-2- Obtain technical simulation parameters from OpenDSS simulation
% Note: here we assume the opertions function can only handle one state at
% a time. In the future this could be vectorized

%Setup storage:
n_states = size(state_list, 1);
backfeed_kWh = NaN(n_states, 1);
net_energy_MWh = NaN(n_states, 1);

%Loop over all states
for s = 1:n_states
    % Extract out state components for easier naming
    pv_peak_fraction = state_list(s, 1);
    storage_kVA = state_list(s, 2);
    storage_kWh = state_list(s, 3);
    
    % Actually run operations
    [backfeed_kWh(s), net_energy_MWh(s)] = SimpleStoragePvOpsRaw(params, pv_peak_fraction, storage_kVA, storage_kWh);
end

%-2- Compute costs 

% Unscaled cost
%  Note: For ISGT problem, baseline technical parameters are for 4 representative days
ops_values = backfeed_kWh/1000 * params.ops_backfeed_cost(t)...
                + net_energy_MWh * params.ops_energy_cost(t);
% Scale to annual
ops_values = ops_values * params.ops_to_annual;

% Finally scale to desired number of years
ops_values = ops_values * ops_scale;

if nargout > 1
    % Here we could run a more detailed simulation if needed, such as for a
    % partially hidden Markov process
    new_full_states = state_list;
end
