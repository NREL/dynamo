function step = ssMcClain(n, old_step, target)
%SSMCCLAIN Recursive Regression McClain Step Size Rule
%
% McClains Formula Provides approximately 1/n convergence at first and then
% levels out at a target value.
%
% Takes a single parameter, the target value.
%
% see Powell (2007) p188 and McClain(1974)
%
% Notes:
%  - It is assumed that the first timestep is n=1 (one-indexed)
%  - When target = 0, convergence is identical to 1/n.
%  - This is a recursive algorithm that requires old_step and does not work
%    for vectors of n

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-11-04 09:40  BryanP      Initial Code
%   2  2011-03-15 10:30  BryanP      Expanded comments & switch to if..else
%   3  2012-01-02 16:20  BryanP      Support for param structure

if isstruct(target)
    target = target.target;
end

if n == 1
    step = 1;
else
    step = old_step ./ (1 + old_step - target);
end
