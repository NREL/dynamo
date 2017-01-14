function [output] = verifyProblemStruct( problem, fieldnames_cell )
%VERIFYDPPROBLEMSTRUCT This function checks the problem for the correct
%structure

if nargin == 1
    fieldnames_cell = {'n_periods', ...
                    'state_set', ...
                    'params', ...
                    'random_items', ...
                    'fTerminalValue', ...
                    'fDecisionSet', ...
                    'fRandomApply', ...
                    'fDecisionApply' };
end
            
            
output = ismember(fieldnames_cell, fieldnames(problem));

if ~all(output)
    missing_fieldnames = fieldnames_cell(~output);
    msg = 'Problem definition is missing some field(s).\n\nThe following appear to be missing:\n';
    msg = [msg, sprintf('\n* %s', missing_fieldnames{:})];
    error('ADP:MissingFieldName', msg)
end

end

