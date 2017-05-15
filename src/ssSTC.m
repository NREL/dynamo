function step = ssSTC(n, old_step, param) %#ok<INUSL>
%SSSTC Recursive Regression search-then-converge Step Size Rule
%
% Search than Converge (STC) allows tuning both the early and late
% convergence properties of the step-size rule. It is effectively a
% combination of harmonic and polynomial rates with an extra parameter to
% provide higher values to early exploration steps
%
% param structure fields are:
%   step2   Base step size aka alpha0 [try 1]
%   a       similar to the harmonic stepsize parameter that decreases
%           convergence rate [try ~10]
%   b       shape parameter that forces low n values to remain close to
%           step2. [try ~1000]
%   beta    polynomial learning rate exponent in range 0.5-1 [try ~1]
%
% Formula is:
%   step(n) = step2 * (b/n + a)/(b/n + a + n^beta)
%
% The formulation used here is a generalized version from presented in
% Powell (2007) Approximate Dynamic Programming, p189
%
% Notes:
%  - It is assumed that the first timestep is n=1 (one-indexed)
%  - vector inputs for n are OK

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   6  2017-05-01 16:12  BryanP      BUGFIX: initialize step properly
%   5  2012-01-02 16:30  BryanP      Always start with step = 1 at n=1
%   4  2011-03-15 10:30  BryanP      Expanded comments
%   3  2011-03-15 10:15  BryanP      Corrected formula to start with n=1 and handle basecase
%   2  2011-03-14 17:40  BryanP      Adapted to handle vectors of n
%   1  2010-11-04 09:40  BryanP      Initial Code

step = ones(size(n));

%Note: we already have step = 1 to handle any n=1

%and any n=2
if param.step2 ~= 1
    valid_idx = n==2;
    step(valid_idx) = param.step2;
    valid_idx = n>2;
else
    valid_idx = n>1;
end

%handle general case
n = n-1;
step(valid_idx) = param.step2 .* (param.b./n(valid_idx) + param.a)./...
             (param.b./n(valid_idx) + param.a + n(valid_idx).^param.beta);
