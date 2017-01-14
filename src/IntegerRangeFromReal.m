function [integer_max, offset, step_size] = IntegerRangeFromReal( min_max, step_size, continuous_chunks )
%INTEGERRANGEFROMREAL Convert a range of real values to integers
%  
% integer_max = IntegerRangeFromReal( min_max, step_size)
%   Returns the number of steps needed to span the space identified by
%   the two row matrix min_max = [min; max] with step_size. 
%
%   step_size=0 is assumed to be continuous and step sizes for these
%   continuous are computed such that the continuous space is divided into
%   chunks (identified by the range 1:11)
%
% ... = IntegerRangeFromReal(min_max, step_size, continuous_chunks)
%   specify the number of chunks to use for continuous values
%
% [integer_max, offset] = IntegerRangeFromReal(...)
%   also return the rounded version of minimum values to use as the offset
%
% [integer_max, offset, step_size] = IntegerRangeFromReal(...)
%   also return the actual step_sizes used (including continuous functions)
%
% The original range can be reconstructed as
%   min_max = zeros(2, length(integer_max));
%   min_max(2,:) = integer_max .* step_size;
%   min_max = bsxfun(@plus, min_max, offset);
%
% Example:
%     >> [integer_max, offset, step_size] = IntegerRangeFromReal([0 1 1.5 3 3.2; 2 2 3 4.2 4], [.5 .5 .3 .2 0])
% 
%     integer_max =
%          4     2     5     6     9
%     offset =
%              0    1.0000    1.5000    3.0000    3.2800
%     step_size =
%         0.5000    0.5000    0.3000    0.2000    0.0800
%
% adapted from a piece of SampleNdRange (v4) by Bryan Palmintier 2016
%
% see also:
%  RealRangeFromInteger

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2016-07-08 00:40  BryanP      Adapted from code in SampleNdRange v4

        if nargin < 3
            continuous_chunks = 10;
        end
        
        % Compute step sizes for any continuous values
        if any(step_size == 0)
            continuous_mask = (step_size == 0);
            step_size(continuous_mask) = ...
                diff(min_max(:, continuous_mask), 1, 1) ./ continuous_chunks;
        end
            
        %Convert the discrete sample range to a set of integers
        d_min = round(min_max(1,:)./step_size, 0);
        d_max = round(min_max(2,:)./step_size, 0);

        %Store range and add one to include both end points.
        integer_max = d_max - d_min;
        offset = d_min .* step_size;
end

