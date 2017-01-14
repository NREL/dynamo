function step = ssSTC2C(n, old_step, param)
%SSSTC2C Recursive Regression search-then-converge-to-constant Step Size Rule
%
% This step-size algorithm, attempts to combine the strong points of:
%   - Harmonic: by slowing the rate of convergence,
%   - McClain: by converging to a constant rather than zero, and
%   - Search than Converge (STC): by allowing a delay before convergence
%
% The resulting formula contains a mix of recursive and deterministic terms
% but behaves similarly to the corresponding STC formula, though with a
% slightly faster initial fall off.
%
% param structure fields are:
%   step2   Base step size aka alpha0 [try 1]
%   a       similar to the harmonic stepsize parameter that decreases
%           convergence rate [try ~1]
%   b       shape parameter that forces low n values to remain close to
%           step2. [try ~1000]
%   target  asymptocially approach this target value [try ~0.002]
%
% Formula is:
%   step(n) = step2 * (b/(n^2) + a) / (b/(n^2) + (a*step2-target)/step(n-1) + 1)
%
% This algorithm was developed by Bryan Palmintier by combining the
% features of search than converge and McClains formulas. For
% further reading see Powell (2007) Approximate Dynamic Programming,
% p187-188
%
% Note: It is assumed that the first timestep is n=1 (one-indexed)
%
% See also ssSTC, ssMcClain, ssHarmonic, ssHarmonic2C, ssSTC2Crecur

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-11-04 09:40  BryanP      Initial Code based on ssSTC v3 & ssMcClain v2
%   2  2012-01-03 15:25  BryanP      Added step2 parameter
%   3  2012-01-02 16:30  BryanP      Always start with step = 1 at n=1

if n == 1
    step = 1;
elseif n == 2 && param.step2 ~= 1
    step = param.step2;
else
    step = param.step2 .* (param.b./(n.^2) + param.a)./...
        (param.b./(n.^2) + (param.a*param.step2 - param.target)./old_step + 1);
end
