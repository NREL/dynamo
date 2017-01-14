function idx_list = mat2ind(matrix_size, subscript_list)
% MAT2IND Convert a list of subscripts to linear indices
%
% Usage: idx_list = mat2ind(matrix_size, subscript_list)
%
% Converts a list of subscripts (one set per row) into linear indexes for
% use with multi-dimensional arrays of unknown or arbitrary numbers of
% dimensions.
%
% Parameters (must be column vectors)
%   matrix_size     The size of the matrix being indexed, as returned by
%                    size(matrix)
%   subscript_list  A matrix list of subscripts with one row per element to
%                    convert to a linear index.
%
% Returns
%   idx_list   A column array of linear indices
%
%
% see also:
%   ind2mat, sub2ind, ind2sub
%
% originally by Bryan Palmintier 2011

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-03-24 15:40  BryanP      Extracted & renamed from faDiscrete helper
%   2  2011-03-29 11:10  BryanP      Fixed to handle lists of states

%%
% Compute the valid number of dimensions so we can ignore any subscripts
% that are beyond this number of dimensions. This can occur when the matrix
% being index has inadvertantly been reduced because trailing dimensions
% have size = 1
valid_dims = length(matrix_size);
if any(subscript_list(:,valid_dims+1:end) > 1)
    warning('Mat2Ind:TooManySubscripts', 'Number of subscripts (>1) exceeds matrix_size')
end

% convert states matrix into a collection of equally sized
% column vectors and save into a row cell vector to be used as
% our list of subscripts
idx_list = mat2cell(subscript_list(:,1:valid_dims), size(subscript_list,1), ones(1,valid_dims));

% now extract the indexes from this collection of subscripts
idx_list = sub2ind(matrix_size,idx_list{:});
end
