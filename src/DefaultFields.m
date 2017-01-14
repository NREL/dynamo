function structure = DefaultFields(structure, field, default)
% SETDEFAULTFIELDS sets/parses missing defaults in an options structure
%
% structure = DefaultFields(structure, field, default) Checks that the
%  fieldname (string) field is contained in structure, if not it adds
%  structure.field and sets to the value default.
%
%  Alternatively, field can be a cell vector of field names and default an
%  equal length cell vector of default values
%
% structure = DefaultFields(structure, default_table) accepts a 2 column
%  cell array, default_table, with field name strings in the first column
%  and default values in the second
%
% Hint:
%  - for lightweight option parsing, use:
%      opt = DefaultFields(struct(varargin{:}),cell_table_of_defaults)
%
% See also: DefaultOpts, inputParser
%
% originally by Bryan Palmintier 2011

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   6  2016-10-06 11:27  BryanP      Renamed to simply DefaultFields to avoid confusion with set classes based on AbstractSet 
%   5  2012-06-25 23:45  BryanP      BUGFIX: fix handle of empty field/default
%   4  2012-06-13 13:05  BryanP      BUGFIX: corrected handling of cell vector inputs
%   3  2012-06-09 22:40  BryanP      Added support cell array options structure
%   2  2012-01-21 10:20  BryanP      Added support for single default_table
%   1  2012-01-13 16:00  BryanP      Extracted from CapPlanDpInit v71

if isempty(structure)
    structure = struct();
elseif iscell(structure)
    %Ensure we have two columns
    if isvector(structure) && length(structure) > 2
        structure = reshape(structure,2,[])';
    end
    %Convert to a structure
    structure = cell2struct(structure(:,2),structure(:,1));
end

if nargin < 3
    if nargin < 2 || isempty(field)
        return
    end
    default = field(:,2);
    field = field(:,1);
end

if not(iscell(field))
    %handle scalar case
    if not(isfield(structure, field))
        structure(1).(field) = default;
    end
else
    %handle vector/cell array case
    for f = 1:length(field)
        %First check if the value exists at all
        if not(isfield(structure, field{f}))
            structure(1).(field{f}) = default{f};
        end
    end
end
end
