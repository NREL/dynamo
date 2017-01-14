function step = ssHarmonic(n, old_step, a)
%SSHARMONIC Recursive Regression harmonic step size
%
% The generalized harmonic step sizes allow for a tunably slower
% relative to 1/n
%
% parameter is a single scaler, a, used to decreas the convergence rate.
% [try ~10]. It can also be passed a a structure in the 'a' field as in
% param.a
%
% Formula is:
%   step(n) = a/(a+n+1)
%
% The formulation is based on eq 6.18 from Powell (2007) Approximate Dynamic
% Programming, p187
%
% Notes:
%  - It is assumed that the first timestep is n=1 (one-indexed)
%  - vector inputs for n are OK

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-03-15 10:00  BryanP      Initial Code
%   2  2012-01-02 16:20  BryanP      Support for param structure

if isstruct(a)
    a = a.a;
end

step = a./(a+n-1);