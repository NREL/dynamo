function [ decision_set, params ] = SimpleStoragePvDecisionSet(params, t, pre_decision_state) 
%SIMPLESTORAGEPVDECISIONSET Produces a set of decisions for the given pre decision state
% Produces a set of decisions for the given pre state (pre_s) at time t
% for use with DP/ADP shared problem structure. 

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   4  2017-04-26 05:35  BryanP      BUGFIX: implement with setCombinWithLimits (finally available)  
%   3  2016-10-27 12:35  BryanP      Added t as second input 
%   2  2016-10-21 15:47  BryanP      Overhaul for new problem structure 
%   1  2010              BryanP      Original version from MIT

%-- Decision options
% First all possible storage combinations
storage_ratings = allcomb(params.storage_kVA{t+1}, params.storage_kWh{t+1});

% Then filter for valid combinations (storage_kWh ./ storage_kVA);
storage_hrs = storage_ratings(:, 2) ./ storage_ratings(:, 1);
ok_idx = and(storage_hrs >= params.storage_hr_limits(1), storage_hrs <= params.storage_hr_limits(2));
% At this point any [0 0] states will be not_ok b/c divide by zero = NaN
% so add back in these [0 0] as valid
ok_idx = or(ok_idx, and(storage_ratings(:, 2) == 0, storage_ratings(:, 1) == 0));
storage_ratings = storage_ratings(ok_idx, :);

% Now compute correspoinding decisions    
decisions = bsxfun(@minus, storage_ratings, pre_decision_state(2:end));

decision_set = setList(decisions);

end

