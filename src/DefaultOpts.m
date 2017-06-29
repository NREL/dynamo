function [options, unused, match_map] = DefaultOpts(opt_in, defaults, varargin)
% DefaultOpts updates defaults based on valid user supplied options
%
% options = DefaultOpts(opt_in, defaults)
%  Checks each field in opt_in against those in defaults. If the field name
%  matches, the value from opt_in is used instead of the default and
%  returned in an updated options structure. If some of the opt_in s are
%  not found in the defaults, a warning is issued. Both opt_in and defaults
%  can be any combination of:
%    1) a structure (scalar)
%    2) a 2 column cell array with field name strings in the first column
%       and default values in the second.
%    3) a cell list of option value pairs, such as contained in varargin
%
% [options, unused] = DefaultOpts(opt_in, defaults)
%  Also returns a cell vector of unused option name strings
%
% [options, unused, match_map] = DefaultOpts(opt_in, defaults)
%  Also returns a boolean map corresponding to which opt_in field names
%  were actually used.
%
% [...] = DefaultOpts(opt_in, defaults, 'warn_unused', false)
%  Suppresses the warning for unmatched fieldnames. In conjuction with
%  match_map, allows the caller to handle field matches across multiple
%  calls to DefaultOpts.
%
% Note: DefaultOpts is very similar to DefaultFields. The difference is
% that DefaultOpts only includes those fields found in defaults and ignores
% any extra fields from opt_in, while DefaultFields keeps all of the
% opt_in fields. As a result, DefaultOpts is more efficient when the
% defaults list is longer than opt_in.
%
% Hint:
%  - for option parsing with unused checks, use:
%      [opt, unused] = DefaultOpts(varargin,cell_table_of_defaults)
%
% See also: DefaultFields, inputParser
%
% originally by Bryan Palmintier 2011

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  10  2017-06-14 20:11  BryanP      Always show unused options warning 
%   9  2017-06-14 06:33  BryanP      BUGFIX add  support for nested cell arrays 
%   8  2017-06-14 06:17  BryanP      BUGFIX when handling single item structs for opt_in 
%   7  2016-10-06 11:27  BryanP      Renamed to simply DefaultOpts to avoid confusion with set classes based on AbstractSet 
%   6  2012-06-27 15:55  BryanP      UPDATE: support intial nested struct + additional pairs in a cell array
%   5  2012-06-26 14:35  BryanP      UPDATE: support for struct nested in a cell array (such as from varargin)
%   4  2012-06-25 21:35  BryanP      RENAMED SetDefaultOpts (was SetOptions) to avoid conflict with built-in (svn: 86)
%   3  2012-06-13 13:05  BryanP      BUGFIX: corrected handling of cell vector inputs
%   2  2012-06-13 10:20  BryanP      Added unused & match_map
%   1  2012-06-09 12:00  BryanP      Adapted from SetDefaultFields v2

my_defaults = { 'warn_unused'   true
              };
my_opt = DefaultFields(varargin, my_defaults);

%Provide a default for early returns
match_map = [];
unused = {};

if isempty(defaults)
    options = [];

% [1] Copy over defaults, converting to structure if needed
elseif iscell(defaults)
    %Check divisable by two to have sufficient name/val pairs
    if mod(numel(defaults), 2) ~= 0
        error('ADP:DefaultOpts:WrongLength', ...
            'defaults must have name-value pairs (numel in %d not divisable by 2)', numel(defaults))
    end
    %Ensure we have two columns
    if isvector(defaults) && length(defaults) > 2
        defaults = reshape(defaults,2,[])';
    end
    %Convert to a structure
    options = cell2struct(defaults(:,2),defaults(:,1));
elseif isstruct(defaults)
    options = defaults;
else
    error('ADP:DefaultOpts:InvalidType','Defaults must be a cell array or structure')
end

% [2] Setup opt_in
if isempty(opt_in)
    % (A) if it is empty, we are done (since we already copied/converted
    % defaults to options
    return
elseif iscell(opt_in)
    %Check divisable by two neded to have sufficient name/val pairs
    if mod(numel(opt_in), 2) ~= 0
        %If not, maybe the first element is a structure we can use
        if isstruct(opt_in{1})
            %In which case use it by merging with defaults. We let the
            %rest of the code handle other option pairs.
            if length(opt_in) > 1
                if nargin < 3
                    varargin = {};
                end
                options = DefaultOpts(opt_in{1}, defaults, varargin{:});
                %Remove our structure so that we now have valid pairs
                opt_in(1) = [];
            else
                %Note: this is also our recursive basecase
                opt_in = opt_in{1};
            end
        %Or maybe it is a nested structure
        elseif iscell(opt_in{1})
            opt_in = opt_in{1};
        else
            %Otherwise we have problems
            error('ADP:DefaultOpts:WrongLength', ...
                'opt_in must have name-value pairs (numel in %d not divisable by 2)', numel(opt_in))
        end
    end
    % (B) convert cell to structure if required
    %If here, have to be divisible by 2 so Ensure we have two columns
    if isvector(opt_in) && length(opt_in) > 2
        opt_in = reshape(opt_in,2,[])';
    end
    %Convert to a structure
    if not(isstruct(opt_in))
        opt_in = cell2struct(opt_in(:,2),opt_in(:,1));
    end
elseif not(isstruct(opt_in))
    error('ADP:DefaultOpts:InvalidType', 'opt_in must be a structure, cell, or empty')
end

% Initialize match map
f_names = fieldnames(opt_in);
n_opt_in = length(f_names);

if isempty(options)
    match_map = false(1, n_opt_in);
    unused = f_names;
else
    match_map = true(1, n_opt_in);

    % [3] Overwrite default values with those from opt_in
    for f = 1: n_opt_in
        f_str = f_names{f};
        %Only copy those with matching fields in defaults
        if isfield(options, f_str)
            options(1).(f_str) = opt_in.(f_str);
        else
            unused{end+1} = f_str; %#ok<AGROW> b/c don't know how big it will be
            match_map(f) = false;
        end
    end

    if isempty(fieldnames(options))
        options = [];
    end
end

if not(isempty(unused)) && my_opt.warn_unused
    unused_str = sprintf('%s ', unused{:});
    unused_str = strtrim(unused_str);
    warning('ADP:DefaultOpts:UnusedOpt','Unused parameters: %s', unused_str)
end

end %Main function
