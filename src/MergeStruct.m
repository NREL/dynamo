function out = MergeStruct(s1, s2)
%Simply merges two structures, with preference for fields defined in later
%listed structure.
%
% Why is this not in MATLAB?

fields = fieldnames(s1);
for f = 1:length(fields);
    out.(fields{f}) = s1.(fields{f});
end

fields = fieldnames(s2);
for f = 1:length(fields);
    out.(fields{f}) = s2.(fields{f});
end