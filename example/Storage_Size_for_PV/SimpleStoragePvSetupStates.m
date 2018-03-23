function states = SimpleStoragePvSetupStates(params)
% SIMPLESTORAGEPVSETUPSTATES Compute list of possible state combinations over time 
%
% states = SimpleStoragePvSetupStates(params_in)
%
% Notes:
%  -- One state per row
%  -- column order: pv_pen storage_kVA storage_kWh
%
% See also: SimpleStoragePvSetupParams, SimpleStoragePvSetupProblem

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   2  2017-07-16 21:24  BryanP      Completed initial parameter, etc. definition 
%   1  2017-07-16 11:19  BryanP      Initial code 

%TODO: duplicate/extend required parameters

states = cell(1, params.n_periods + 1);
states{1} = [params.pv_state_set{1}(1), params.storage_kVA{1}(1), params.storage_kWh{1}(1)];
for t = 2:params.n_periods + 1 
    % First all possible combinations
    states{t} = allcomb(params.pv_state_set{t}, params.storage_kVA{t}, params.storage_kWh{t});

    % Then filter for valid combinations (storage_kWh ./ storage_kVA);
    storage_hrs = states{t}(:, 3) ./ states{t}(:, 2);
    ok_idx = (storage_hrs >= params.storage_hr_limits(1)) & (storage_hrs <= params.storage_hr_limits(2));
    states{t} = states{t}(ok_idx, :);
end 