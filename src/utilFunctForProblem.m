function f_handle = FunctForProblem(problem, f_string, default)
% FUNCTFROMPROB Establish function (pointer) for a problem structure/object
%
% f_handle = FunctFromProb(problem, f_string, default)
%
% Helper function to setup required function handles based on either
% structure or object based problems, with a default if not specified by
% the problem

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2012-04-16 17:35  BryanP      Extracted from parTD1 v3
%   1  2012-04-16 20:15  BryanP      Throw error when default not allowed

if isstruct(problem) && isfield(problem, f_string) && not(isempty(problem.(f_string)))
    f_handle = problem.(f_string);
elseif isobject(problem) && ismethod(problem, f_string)
    f_handle = @(varargin)problem.(f_string)(varargin{:});
elseif isnumeric(default) && not(isempty(default)) && isnan(default)
    error('ADP:FunctForProblem:MustDef', 'Problem must define %s', f_string)
else
    f_handle = default;
end

end
