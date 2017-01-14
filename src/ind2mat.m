function subscript_list = ind2mat(matrix_size, idx_list)
% IND2MAT Convert a list of linear indices to subscripts as a matrix
%
% Usage: subscript_list = ind2mat(matrix_size, idx_list)
%
% Converts a vector of linear indexes into a list of subscripts (one set
% per row) for use with multi-dimensional arrays of unknown or arbitrary
% numbers of dimensions.
%
% Parameters (must be column vectors)
%   matrix_size     The size of the matrix being indexed, as returned by
%                    size(matrix)
%   idx_list   A column array of linear indices
%
% Returns
%   subscript_list  A matrix list of subscripts with one row per element to
%                    convert to a linear index.
%
% see also:
%   mat2ind, sub2ind, ind2sub
%
% originally by Bryan Palmintier 2011

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-03-25 19:45  BryanP      Extracted & renamed from faDiscrete

% preconstruct a cell array to receive the multiple
% dimensions of our subscript
subscript_list = num2cell(zeros(1,length(matrix_size)));

% Convert the linear find index to subscripts and assign to
% each element of the cell array
[subscript_list{:}] = ind2sub(matrix_size, reshape(idx_list,[],1));

subscript_list = horzcat(subscript_list{:});
end
