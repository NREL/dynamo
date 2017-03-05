function f_handle = utilFunctForProblem(problem, f_string, default)
% UTILFUNCTFROMPROB Establish function (pointer) for a problem structure/object
%
% f_handle = utilFunctFromProb(problem, f_string, default)
%
% Helper function to setup required function handles based on either
% structure or object based problems, with a default if not specified by
% the problem
% 
% For use with the adp* family of functions.

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   3  2017-03-04 23:35  BryanP      Prepend "util" to function name
%   2  2012-04-16 20:15  BryanP      Throw error when default not allowed
%   1  2012-04-16 17:35  BryanP      Extracted from parTD1 v3

if isstruct(problem) && isfield(problem, f_string) && not(isempty(problem.(f_string)))
    f_handle = problem.(f_string);
elseif isobject(problem) && ismethod(problem, f_string)
    f_handle = @(varargin)problem.(f_string)(varargin{:});
elseif isnumeric(default) && not(isempty(default)) && isnan(default)
    error('ADP:utilFunctForProblem:MustDef', 'Problem must define %s', f_string)
else
    f_handle = default;
end

end
