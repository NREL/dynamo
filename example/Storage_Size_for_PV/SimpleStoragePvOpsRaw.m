function [backfeed_kWh, net_energy_MWh] = SimpleStoragePvOpsRaw(params, pv_peak_fraction, storage_kVA, storage_kWh)
% SIMPLESTORAGEPVOPSRAW Raw Operations costs for simple storage + PV problem
%
%   Handles a single state. Use params to pass options that may change with
%   various simulations. For example, the feeder name and list of days to
%   simulate.

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-07-17 08:43  BryanP      Adapted MultiInvOps v3

% Psuedo-code:
% Handle memoization (?)

ipaddress = 'tcp://localhost:5556';

if params.run_opendss
    for d_idx = 1: length(params.days_to_simulate)
        date_str = params.days_to_simulate{d_idx};
        j.pv_peak_fraction = pv_peak_fraction;
        j.storage_kVA = storage_kVA;
        j.storage_kWh = storage_kWh;
        j.date_str = date_str;
        data = savejson('', j, 'Compact', 1);
        results = zmq_client(data, ipaddress);
        backfeed_kWh = results.backfeed_kWh;
        net_energy_MWh = results.net_energy_MWh;
    end
    % Add up across days keeping by hour resolution

else
    % Use lookup table
    [~, state_idx] = ismembertol([pv_peak_fraction, storage_kVA, storage_kWh], ...
                            params.stored_ops(:,1:3), 'ByRow', true);
    %Keep first if multiple matches
    state_idx = state_idx(1);
    backfeed_kWh = params.stored_ops(state_idx, 4);
    %Sum net energy in kWh and convert to MWh
    net_energy_MWh = sum(params.stored_ops(state_idx, 5:end))/1000;
end

% Sum net_energy_WWh across hours

