function [data, table_specifier] = ReadGamsTable(filename, has_header, delim, ...
    data_in_rows, data)
%READGAMSTABLE Read data from GAMS input file (either gms or delimited)
%
% [data, table_specifier] = ReadGamsTable(filename, has_header, delim, ...
%    data_in_rows, data)
%
% Reads the first table into a struct array with field names taken from
% the table labels
%
% For data_in_rows = true:
%  - The array of structures allows us to easy access to both all of a
%    generators parameters (by using data.list(INDEX) or to access all of
%    the same attribute across generators using data.list(:).ATTRIBUTE
%  - But, it does not let us easily look up which index corresponds to
%    which name. Since GAMS seems to like to rearrange the order of the
%    names, we will also create a mapping structure to lookup by name, such
%    as data.map.NAME which will return the corresponding index for
%    use with data.list(INDEX)

%
% Note: supports GAMS comments (* in first column) for all input data
%
% originally by Bryan Palmintier 2011

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-01-18 15:47  BryanP      Original draft version
%   2  2011-01-20 09:45  BryanP      Fully working version!
%   3  2011-01-20 20:45  BryanP      Store filename in data_file field
%   4  2011-01-21 21:05  BryanP      Added support for data_in_rows
%   5  2011-01-21 22:00  BryanP      Able to append to existing data struct
%   6  2011-01-22 24:10  BryanP      Store file name and path separately
%   7  2011-11-07 17:25  BryanP      New default delim, any # of: spaces, tabs, and commas
%   8  2011-11-07 17:25  BryanP      Fixed errors related to first column & default delim

%Setup flags based on inputs
has_header = nargin < 2 || isempty(has_header) || has_header;

use_default_delim = nargin < 3 || isempty(delim);
if use_default_delim
    delim = ', \t';
end

data_in_rows = not(nargin <4 || isempty(data_in_rows)) && data_in_rows;

if nargin <5 || not(isstruct(data))
    data = struct([]);
end

% 0) Open the file and suck it into a giant string
f_as_string = fileread(filename);

%Also store the file name & path for later
[pathstr, name, ext] = fileparts(filename);
%Have to include (1) if working with a blank structure
data(1).data_path = pathstr;
%Note: the '.' is included in ext
data.data_file = [name ext];


% 1) Remove the "header" = to all into non-table preamble, skipping
%   everything upto the table keyword

if has_header
    %Build up regular expression search pattern:
    % a) look for the keyword table (table) at the beginning of a line (ie
    % after a newline, \n) possibly preceed by whitespace (\s*) and
    % definitely followed by one or more whitespace characters (\s+)
    search_pattern = '\n\s*table\s+';
    % b) match the table specifier consisting of the table name and
    % associated set definition. The surrounding ()'s specify that we want
    % to save this output. ?<NAME> defines the field name to which to save
    % it. \w specifies the name must start with a word character.
    % [\w\d,\(\s]+ specifies that any word, number, or space character,
    % plus a comma and the open parentheses make up the middle part of the
    % table specifier. \) tells us to end the specifier with a close paren.
    % We have to escape the parentheses with the backslash here (and
    % earlier) because they have special meaning in regular expressions.
    search_pattern = [search_pattern '(?<specifier>\w[\w\d,\(\s]+\))'];
    % c) Now skip to the end of the line, by matching any number of non-newline
    % characters ([^\n]* and finally a newline \n
    search_pattern = [search_pattern '[^\n]*\n'];
    % d) Now extract everything up until the first semicolon as the data
    % table itself. As in b), the (?<NAME>)s save what is inside to the
    % named field. [^;]* matches any number of non-semicolons. Then we
    % check for (but don't save) the final semicolon).
    %
    % Implementation Note: this assumes no semi-colons in the comments.
    search_pattern = [search_pattern '(?<data>[^;]*);'];

    %Finally, we apply our search pattern to the file and extract only the
    %first data table to struct array with our specified field names.
    table = regexpi(f_as_string, search_pattern, 'once','names');
else
    table.data = f_as_string;
    table.specifier = [];
end

% 2) Remove any comment lines (use the annoying GAMS convention of starting
% comments with a * at the start of the line and then also using * for
% ranges of field names (valid data) and multiplication and... ergh
table.data = regexprep(table.data, '\n\*[^\n]*', '');

% 3) Now Identify the field names in the first row

% Extract the first row
[labels, data_start] = regexp(table.data, '([^\n]+\n)', 'once', 'tokens','end');

if isempty(labels)
    warning('ReadGamsTable:no_labels', 'Can''t find label row. Make sure file has linefeed characters \n at end of each line')
end

% And remove it from the table
table.data(1:data_start)=[];

%See if we specified an alternate name for the first column
first_col = regexp(labels, ['^([^' delim ']+)'], 'once', 'tokens');
%if so, our search string will find it so clear the default, otherwise
%set it to our default
if isempty(first_col) || isempty(first_col{1})
    first_col = 'id';
else
    first_col = [];
end

%And build up our search string for the remaining column labelss based
%around our delimiter
search_pattern = ['\s*([^\r\n' delim ']+)\s*[\n' delim ']'];

fields = regexp(labels, search_pattern,'tokens');
fields = vertcat(first_col,fields{:}{:});

% 4) Read in the first line of data as strings to deterine the column data
% types
% Extract the first row (but leave it in the table)
first_data = regexp(table.data, '([^\n]+\n)', 'once', 'tokens');

first_data = regexp(first_data, search_pattern,'tokens');
first_data = vertcat(first_data{:}{:});

if length(first_data) ~= length(fields)
    warning('ReadGamsTable:HeadDataLenMismatch', 'The labels row and first data row are of different lengths')
end

format_str = '';

for f = 1:length(first_data)
    if isempty(str2num(first_data{f})) %#ok<ST2NM> %use str2num b/c we want to allow NaN as valid input
        format_str = [format_str '%s']; %#ok<AGROW>
    else
        format_str = [format_str '%f']; %#ok<AGROW>
    end
end

% 5) Use textscan to read in the rest of the data using the formats
% determined above
if use_default_delim
    data_table = textscan(table.data,format_str, 'Delimiter', delim, 'MultipleDelimsAsOne', true);
else
    data_table = textscan(table.data,format_str, 'Delimiter', delim);
end

% 6) Convert the results to a structure array
if data_in_rows
    %Add a list and map substructure each with one field for each row
    for r = 1:length(data_table{1})
        field_name = data_table{1}{r};
        data.map.(field_name) = r;
        for f = 1:length(fields)
            val = data_table{f}(r);
            if isnumeric(val)
                data.list(r).(fields{f}) = data_table{f}(r);
            else
                data.list(r).(fields{f}) = data_table{f}{r};
            end
        end
    end
    data.n = r;
else
    %Otherwise, convert to a 1-D structure with a vector for each
    %field/column
    for f = 1:length(fields)
        data.(fields{f}) = data_table{f};
    end
    data.n = length(data_table{1});
end



if nargout > 1
    table_specifier = table.specifier;
end
