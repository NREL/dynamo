function step = ssSTC2Crecur(n, old_step, param)
%SSSTC2Crecur Fully Recursive search-then-converge-to-constant Step Size Rule
%
% This step-size algorithm is a modified version of STC2C that does not
% rely on the iteration number, n, providing better results when resetting
% step sizes in the middle of analysis such as in multi-fidelity. Like
% STC2C it attempts to combine the strong points of:
%   - Harmonic: by slowing the rate of convergence,
%   - McClain: by converging to a constant rather than zero, and
%   - Search than Converge (STC): by allowing a delay before convergence
%
% The resulting formula contains purely recursive terms and behaves
% similarly to the corresponding STC and STC2C formulas. However, this form
% tends to drop even more slowly than the non-recursive STC2C requiring a
% further reduction of b to achieve similar trends to STC (or the basic
% STC2C)
%
% param structure fields are:
%   step2   Base step size aka alpha0 [try 1]
%   a       similar to the harmonic stepsize parameter that decreases
%           convergence rate [try ~1]
%   b       shape parameter that forces low n values to remain close to
%           step2. [try ~10]
%   target  asymptocially approach this target value [try ~0.002]
%
% To compute the formula, we first estimate the current n values using the
% results of backsolving the basic STC formula for n using the quadratic
% formula and then adding one. In this calculation the term (1-1/step)
% occurs frequently. In this case we adjust this value not just by step2
% (as in STCrecur), but also by the target value. These result leading to
% the series of equations:
%   s_fact = 1-1/((step(n-1)-target)/(step2-target))
%   n_est = 1 + (- a*s_fact + sqrt(s_fact^2*a^2 - 4*s_fact*b) )/2;
% and finally
%   step(n) = step2 * (b/n_est^2 + a)/(b/n_est^2 + (a*step2 - target)/step(n-1) + 1)
%
% This algorithm was developed by Bryan Palmintier by combining the
% features of search than converge and McClains formulas. For
% further reading see Powell (2007) Approximate Dynamic Programming,
% p187-188
%
% Note: It is assumed that the first timestep is n=1 (one-indexed)
%
% See also ss1overN, ss1overNrecur, ssMcClain, ssHarmonic, ssHarmonicRecur,
%   ssHarmonic2C, ssSTC, ssSTC2C, ssSTCrecur

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2012-01-02 12:30  BryanP      Initial Code based on ssSTCrecur v2 and ssSTC2C v2
%   2  2012-01-02 16:30  BryanP      Always start with step = 1 at n=1

if n == 1
    step = 1;
elseif n == 2 && param.step2 ~= 1
    step = param.step2;
else
    s_fact = 1 - 1/((old_step-param.target)/(param.step2-param.target));
    n_est = 1 + (- param.a*s_fact + sqrt(s_fact^2*param.a^2 - 4*s_fact*param.b) )/2;

    step = param.step2 .* (param.b./(n_est.^2) + param.a)./...
        (param.b./(n_est.^2) + (param.a*param.step2 - param.target)./old_step + 1);
end
