function t_or_f = utilEqualWithTol(a, b, tol)
% UTILEQUALWITHTOL equality test with tolerance
%
% t_or_f = utilEqualWithTol(a, b, tol)
%
% Helper function for simple relative value compare with divide by zero
% work around
% 
% For use with the adp* family of functions.

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2012-04-16 17:35  BryanP      Extracted from adpTD1 v21

    t_or_f =  abs(a-b)./max(max(max(abs(a),[],2),max(abs(b),[],2)),1e-3) < tol;
end
