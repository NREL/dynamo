function step = ss1overNrecur(n, old_step, param) %#ok<INUSD>
%SS1overNrecur Recursive Regression version of 1/n Step Size Rule
%
% This step-size algorithm is a modified version of 1/n that does not
% rely on the iteration number, n. This provides better results when
% resetting step sizes in the middle of analysis such as in multi-fidelity.
%
% The resulting formula contains purely recursive terms and behaves
% identically to the 1overN formula.
%
% Notes:
%  - It is assumed that the first timestep is n=1 (one-indexed)
%  - This is a recursive algorithm that requires old_step and does not work
%    for vectors of n
%
% See also ss1overN, ssMcClain, ssHarmonic, ssHarmonicRecur,
%   ssHarmonic2C, ssSTC, ssSTC2C, ssSTCrecur, ssSTC2Crecur

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-11-04 09:40  BryanP      adapted from ssMcClain v3

if n == 1
    step = 1;
else
    step = old_step ./ (1 + old_step);
end
