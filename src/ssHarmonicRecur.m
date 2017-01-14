function step = ssHarmonicRecur(n, old_step, a)
%SSHARMONICRecur Recursive Regression harmonic step size
%
% This step-size algorithm is a modified version of Harmonic step sizes
% that does not rely on the iteration number, n. This provides better
% results when resetting step sizes in the middle of analysis such as in
% multi-fidelity.
%
% The resulting formula contains purely recursive terms and behaves
% similarly to the corresponding Harmonic formula by allows a tunably
% slower convergence relative to 1/n.
%
% parameter is a single scaler, a, used to decreas the convergence rate.
% [try ~10]. It can also be passed a a structure in the 'a' field as in
% param.a
%
% This algorithm was developed by Bryan Palmintier by simplifying the
% algorithm in ssHarmonic2C. For further reading see Powell (2007)
% Approximate Dynamic Programming, p187-188
%
% Notes:
%  - It is assumed that the first timestep is n=1 (one-indexed)
%  - This is a recursive algorithm that requires old_step and does not work
%    for vectors of n
%
% See also ssSTC, ssMcClain, ssHarmonic, ssSTC2C, ssHarmonic2C


% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2012-01-02 16:10  BryanP      Adapted from ssHarmonic2C v1

if n == 1
    step = 1;
    return
end

if isstruct(a)
    a = a.a;
end

step = a * old_step ./ (a + old_step);
