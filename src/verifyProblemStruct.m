function [problem] = verifyProblemStruct( problem, required_fields )
%VERIFYDPPROBLEMSTRUCT checks dynamo problem for  correct structure 
%
% Also ensures that all optional fields are included to allow checking for
% them with not(isempty())

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   2  2017-09-24 21:39  BryanP      Cleaned up, added optional field creation as [] 
% many                               Early versions


%% Part 1: error if required fields are not included
if nargin == 1
    required_fields = {'n_periods', ...
                    'state_set', ...
                    'params', ...
                    'random_items', ...
                    'fTerminalValue', ...
                    'fDecisionSet', ...
                    'fRandomApply', ...
                    'fDecisionApply' };
end
            
required_field_found_map = isfield(problem, required_fields);

if  not(all(required_field_found_map))
    missing_fieldnames = required_fields(~required_field_found_map);
    msg = 'Problem definition is missing some required field(s).\n\nThe following appear to be missing:\n';
    msg = [msg, sprintf('\n* %s', missing_fieldnames{:})];
    error('ADP:MissingFieldName', msg)
end

%% Part 2: Set any missing optional fields to []

optional_fields = {'fOpsBeforeDecision'
                'fOpsAfterDecision'
                'fOpsAfterRandom'
                'fCompareDecisionPolicy'
                'fCompareValue'
                'fMapState2Vfun'
                'fOptimalDecision'
                'fRandomSample'
                'fRandomJoint'
                };

for f_idx = 1:length(optional_fields)
    opt_field = optional_fields{f_idx};
    if not(isfield(problem, opt_field))
        problem.(opt_field) = [];
    end
end


end

