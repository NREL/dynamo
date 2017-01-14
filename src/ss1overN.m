function step = ss1overN(n, old_step, params) %#ok<INUSD>
%SS1OVERN simple 1/n Recursive Regression step size
%
% One of the simplest step-size rules: 1/n
%
% Formula is:
%   step(n) = 1/n
%
% The formulation is based on eq 6.18 from Powell (2007) Approximate Dynamic
% Programming, p187
%
% Notes:
%  - It is assumed that the first timestep is n=1 (one-indexed)
%  - vector inputs for n are OK
%
% See also ss1overNrecur, ssMcClain, ssHarmonic, ssHarmonicRecur,
%   ssHarmonic2C, ssSTC, ssSTC2C, ssSTCrecur, ssSTC2Crecur

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-03-15 10:00  BryanP      Initial Code

step = 1./n;