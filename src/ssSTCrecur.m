function step = ssSTCrecur(n, old_step, param)
%SSSTCrecur Fully Recursive search-then-converge Step Size Rule
%
% This step-size algorithm is a modified version of STC2C that does not
% rely on the iteration number, n. This provides better results when
% resetting step sizes in the middle of analysis such as in multi-fidelity.
%
% The resulting formula contains purely recursive terms and behaves
% similarly to the corresponding STC formula.
%
% param structure fields are:
%   step2   Base step size aka alpha0 [try 1]
%   a       similar to the harmonic stepsize parameter that decreases
%           convergence rate [try ~1]
%   b       shape parameter that forces low n values to remain close to
%           step2. [try ~10]
%
% To compute the formula, we first estimate the current n value using the
% results of backsolving the basic STC formula for n using the quadratic
% formula and then adding one. In this calculation the term (1-1/step)
% occurs requently leading to the series of equations:
%   s_fact = 1-1/(step(n-1)/step2)
%   n_est = 1 + (- a*s_fact + sqrt(s_fact^2*a^2 - 4*s_fact*b) )/2;
% and finally
%   step(n) = step2 * (b/n_est + a)/(b/n_est + a + n_est)
%
% This algorithm was developed by Bryan Palmintier after some algebra. For
% further reading see Powell (2007) Approximate Dynamic Programming,
% p187-188
%
% Note: It is assumed that the first timestep is n=1 (one-indexed)
%
% See also ssSTC, ssMcClain, ssHarmonic, ssHarmonic2C

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2012-01-02 12:30  BryanP      Initial Code based on ssSTC2C v1
%   2  2012-01-03 15:25  BryanP      Added step2 parameter
%   3  2012-01-02 16:30  BryanP      Always start with step = 1 at n=1

if n == 1
    step = 1;
elseif n == 2 && param.step2 ~= 1
    step = param.step2;
else
    s_fact = 1 - 1/(old_step/param.step2);
    n_est = 1 + (- param.a*s_fact + sqrt(s_fact^2*param.a^2 - 4*s_fact*param.b) )/2;

    step = param.step2 * (param.b ./ n_est + param.a)./...
        (param.b ./n_est + param.a + n_est);
end
