function [min_max, real_values] = RealRangeFromInteger(integer_max, offset, step_size)
%REALRANGEFROMREAL Identify real valued min_max from integer, step_size, and offset
%
% min_max = RealRangeFromInteger(integer_max, offset, step_size)
%   Returns a two row column of the form [min;max] for corresponding real
%   values for the specified number of steps of step_size starting at
%   offset
%
% [min_max, real_values] = RealRangeFromInteger(...)
%   Also return a cell array with the lists of corresponding values
%
% originallly by Bryan Palmintier 2016
%
% see also:
%  IntegerRangeFromReal

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   2  2017-07-16 17:17  BryanP      Updated for 1-based indexing in integer_max
%   1  2016-07-08 01:30  BryanP      Initial code

% Initialize storage
min_max = zeros(2, length(integer_max));
% Setup 0:Integer range
min_max(2,:) = (integer_max - 1) .* step_size;
% Now adjust offset
min_max = bsxfun(@plus, min_max, offset);

% If desired return cell array with corresponding row values for each step
if nargout > 1
    real_values = cell(1, length(integer_max));
    for dim = 1:length(integer_max)
        real_values{dim} = (0:(integer_max(dim)-1))' .* step_size(dim) + offset(dim);
    end
end