function step = ssHarmonic2C(n, old_step, param)
%SSHARMONIC2C Recursive Regression harmonic step size to constant
%
% This variant on the generalized harmonic step sizes allow for both a
% tunably slower convergence relative to 1/n and convergence to a specifed,
% typically non-zero, step-size. (Similar to McClain)
%
% The param structure must contain the following fields:
%   a       used to decreas the convergence rate. [try ~10]
%   target  asymptotic value to approach as n -> Inf
%
% This algorithm was developed by Bryan Palmintier by combining the
% features of generalized harmonic step-sizes and McClains formula. For
% further reading see Powell (2007) Approximate Dynamic Programming,
% p187-188
%
% Notes:
%  - It is assumed that the first timestep is n=1 (one-indexed)
%  - When target = 0, convergence is identical to ssHarmonic.
%  - This is a recursive algorithm that requires old_step and does not work
%    for vectors of n
%
% See also ssSTC, ssMcClain, ssHarmonic, ssSTC2C


% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-03-15 10:30  BryanP      Adapted from ssMcClain v2

if n == 1
    step = 1;
    return
end
step = param.a * old_step ./ (param.a + old_step - param.target);
