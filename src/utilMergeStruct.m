function out = utilMergeStruct(s1, s2)
% UTILMERGESTRUCT Simply merges two structures, with preference for fields defined in later
% listed structure.
%
%  out = utilMergeStruct(s1, s2)
%    Merges structures s1 and s2. duplicate fields use those from s2
%
% Why is this not in MATLAB?

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   2  2017-06-14 05:35  BryanP      add util prefix and expand comments 
%   1  2012?             BryanP      Initial version 

fields = fieldnames(s1);
for f = 1:length(fields)
    out.(fields{f}) = s1.(fields{f});
end

fields = fieldnames(s2);
for f = 1:length(fields)
    out.(fields{f}) = s2.(fields{f});
end