function crf = CapitalRecoveryFactor(rate, life)
%CapitalRecoveryFactor compute capital recovery factor (annualized capital)
%
% Usage:
%    crf = CapitalRecoveryFactor(rate, life)
%
% originally by Bryan Palmintier 2011

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-06-11 00:05  BryanP      Adapted from CapPlanDpOpsModel v7

    crf = rate .* ((1 + rate).^life)/((1 + rate).^life - 1);
end
