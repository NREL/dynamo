function partit = RandPartition(space_per_item, total_space, varargin)
% RANDPARTITION randomly divides a given space among multiple items of arbitrary size
%
% Returns the partition, not the number of items. Uses efficient symmetric
% Dirichlet with alpha=1 sampling to sample over the the hyperdiagonal (n-1
% simplex) 
%
% Usage
%  partit = RandPartition(space_per_item, total_space, option_pairs)
%
% Notes:
%  -- space_per_item is row oriented. 
%  -- Generates multiple sample partitions by specifying a vector for
%  total_space or a multi-row space_per_item. If both the number of rows must match
%  -- round() used as the default fFinalRound, which means the partit sum
%      maybe a bit higher or lower than the total_space
%  -- Use a space_per_item of 0 for items capable of continuous space assignment.
%
% Option pairs & defaults:
%     'fRound'           @round   % function used for final rounding
%     'lumpy'            false    % encourage more zeros in partition (default=approx. uniform)
%     'lump_factor'      []       % adjusts lumpy distributions. (see below)
%     'lump_max_nz_fract' []       % randomly add zeros to the (un-normalized)
%                                 % sample space to ensure no more than
%                                 % this fraction (often much less) of the
%                                 % result is non-zero
%     'allow_extend'     false    % enables more than just the final item to extend
%
%     %Options for external random number input. Enables use of sobol, etc.
%     % Note: In limited testing Sobol performed (much) worse than pure random
%     'sample_set'       []       % set of (0,1] random values to use for partioning (n_samples x n_items)
%     'round_idx_sample' []       % list of (0,1] random values to use in selecting the element to round
%
% Algorithms:
%  lumpy=false: randomly samples the hyperdiagonal (n simplex). This
%               results in a fairly uniform distribution of partitions.
%  lumpy=true, lump_factor=#: The lump_factor is multiplied times the
%               samples used in the (pre-normalized) basic partition. This
%               increases/decreases the quantity of each item, thereby
%               filling up the total space faster/slower to introduce
%               additional lumpiness. This results in a highly non-linear
%               relationship between the lump_factor and number of lumps.
%               Some reference values for lump_factor are:
%                   1       same as non-lumpy
%                   2       about n_dims/2 lumps
%                   n_dims  about 2 lumps
%  lumpy=true, lump_max_nz_fract=#: The non-normalized


% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2012-03-26        BryanP      Created in parallel with SampleCPDpState
%   2  2012-06-21 12:12  BryanP      BUGFIX: Prevent negative space
%   3  2012-06-21 20:42  BryanP      FEATURE: add lumpy & lump_factor support
%   4  2012-06-25 21:35  BryanP      updated for renamed SetDefaultOpts (was SetOptions in adp svn <86)
%   5  2012-06-28        BryanP      Support for single item
%   6  2012-06-29        BryanP      Additional tuning: lump_rand_nz, allow_extend
%   7  2012-07-02 23:50  BryanP      MAJOR overhaul:
%                                       -- use expon sampling (not sort)
%                                       -- streamlined lump creation
%                                       -- enable external random values
%                                       -- additional vectorization
%   8  2012-07-03 11:35  BryanP      BUGFIX: need full random sort for rounding

% [1] Setup
% optional parameters and their defaults
defaults = {
            'fRound'           @round   % function used for final rounding
            'lumpy'            false    % encourage more zeros in partition (default=approx. uniform)
            'lump_factor'      []       % pre-multiply random numbers by this
            'lump_max_nz_fract' []       % randomly add zeros to the (un-normalized)
                                        % sample space to ensure no more than
                                        % this fraction (often much less) of the
                                        % result is non-zero
            'allow_extend'     false    % enables more than just the final item to extend

            %Options for external random number input. Enables use of sobol, etc.
            % Note: In limited testing Sobol performed (much) worse than pure random
            'sample_set'       []       % set of (0,1] random values to use for partioning (n_samples x n_items)
            };
opt = DefaultOpts(varargin, defaults);

n_items = size(space_per_item, 2);
n_samples = max(size(space_per_item, 1), size(total_space,1));

sort_for_round = any(space_per_item > 0);

%duplicate space per item if needed
if size(space_per_item,1) == 1
    space_per_item = repmat(space_per_item, n_samples, 1);
end


if n_items > 1
    % [2] Uniform sample over the hyperdiagonal (n-1 simplex) using the symmetric
    % Dirichlet with alpha=1. Random samples of this can be genereated using
    % normalized exponentials.

    %  (A) create exponential samples
    if isempty(opt.sample_set)
        partit = exprnd(1,n_samples, n_items);
    else
        assert(isequal([n_samples, n_items], size(opt.sample_set)), ...
            'ADP:RandPartition:SizeMismatch', ...
            'Optional external sample_set must match output size (%d x %d)', ...
                n_samples, n_items)
        partit = expinv(opt.sample_set,1);
    end

    %  (B) zero out required fraction
    if opt.lumpy && not(isempty(opt.lump_max_nz_fract))
        zero_out_mask = rand(numel(partit),1) < opt.lump_max_nz_fract;
        partit(zero_out_mask) = 0;
    end

    %  (C) Normalize samples
    partit = bsxfun(@rdivide, partit, sum(partit, 2));

    %  (D) Adjust normalized samples for lumpiness if desired
    if opt.lumpy && not(isempty(opt.lump_factor))
        partit = partit * opt.lump_factor;
    end


    %  (E) scale up from [0 1] to fill the entire space
    %
    % Note: bsxfun allows multiplying column by column without looping or
    % creating a redundant repmat()
    partit = bsxfun(@times, partit, total_space);

    % [3] Sort for rounding
    if sort_for_round
        % (B) Create re-order map
        %  ii> Create list (by row) of indicies that move the desired item
        % to the end for rounding
        new_order = cell2mat(arrayfun(@randperm, repmat(n_items,n_samples,1), ...
                                'UniformOutput', false));

        %  iii> convert to linear indices (matlab goes column by column)
        new_order = n_samples * (new_order - 1);
        new_order = bsxfun(@plus, new_order, (1:n_samples)');

        % (C) actually reorder space per item & weights
        space_per_item = space_per_item(new_order);
        partit = partit(new_order);
    end

    % [4] Fill out the space
    space_filled = zeros(n_samples, 1);
    for item = 1:n_items
        % limit this item to the remaining space if needed
        if item == n_items
            partit(:,item) = max(0,total_space - space_filled);
        else
            partit(:,item) = max(0,min(total_space - space_filled, partit(:,item)));
        end

        %round to item size where needed, if space=0 (continous) then there is
        %nothing more to do
        to_round = (space_per_item(:,item) > 0);

        %if we allow extension, then we don't update the remaining space based
        %on the ammount used by this item.
        if opt.allow_extend
            space_filled = space_filled + partit(:,item);
        end

        if item < n_items
            partit(to_round, item) = RoundTo(partit(to_round, item), space_per_item(to_round, item));
        else
            partit(to_round, item) = opt.fRound(partit(to_round, item)./space_per_item(to_round, item))...
                                        .* space_per_item(to_round, item);
        end

        % To stay within the limits, update the currently used space
        if not(opt.allow_extend)
            space_filled = space_filled + partit(:,item);
        end
    end

    % [5] clean up
    if sort_for_round
        %return to original order
        partit(new_order)=partit;
    end

else
    % [6] Simple case for single item
    partit = total_space;
    to_round = (space_per_item > 0);
    partit = opt.fRound(partit(to_round)./space_per_item(to_round))...
                                    .* space_per_item(to_round);
end
