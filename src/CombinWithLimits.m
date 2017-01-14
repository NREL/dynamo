function [item_count_list, value_list, sum_list, opt] = CombinWithLimits(space_per_item, ...
                               total_space, item_max, min_space_used, ...
                               item_min, fRoundMax, fRoundMin)
% CombinWithLimits list space constrained combinations for multi-product inventory
%
% Usage: [item_count_list, value_list, sum_list] = CombinWithLimits(space_per_item, ...
%                               total_space, item_max, min_space_used, ...
%                               item_min, fRoundMax, fRoundMin)
%
% Recursively generates a list of possible space constrained item
% combinations to be used for state space in the multi product inventory
% and related problems.
%
% Parameters (must be column vectors)
%   space_per_item    amount of storage space required for each item. Must
%                      be a column vector with one row for each item
%   total_space       total storage space available in warehouse (scalar)
%   item_max          (optional) vector of per_item max number of items
%                      default: it is OK to fill the entire space with a
%                      single item. i.e. fRoundMax(total_space./space_per_item)
%   min_space_used    (optional) minimum total inventory in warehouse.
%                      default: 0
%   item_min          (optional) minimum inventory for each item (vector)
%                      default: zero for each item
%   fRoundMax         (optional) function handle to use for rounding
%                      non-integer combinations in comparision to the
%                      maximum space defined by total_space. default:
%                      @floor. If you consider boxes for each item in a
%                      spherical warehouse, the use of floor ensures that
%                      all boxes fit within the sphere, while @ceil would
%                      allow the corners of boxes to extend beyond the
%                      sphere, thereby maintaining a spherical shape, just
%                      larger than the reference sphere.
%   fRoundMin         (optional) function handle to use for rounding
%                      non-integer combinations when establishing the
%                      minimum space as specified by min_space_used.
%                      default: complements fRoundMax as defined below:
%                         fRoundMin   fRoundMax    Interpretation
%                           ceil        floor       Strictly within limits (default)
%                           round       round       At least half within limits
%                           floor       ceil        Any portion within limits
%
% Returns
%   item_count_list   2-D matrix with num_itms columns. each row
%                     corresponds to a valid combination of inventory
%                     levels with columns containing the corresponding
%                     inventory level for each product type
%   value_list        2-D matrix with the same layout as item_count_list
%                     that contains the total space per item rather than
%                     count per item
%   sum_list          column vector with the row sum of value_list
%
% Notes:
%  - This function has uses far beyond the multi-product inventory problem.
%     Ex- Electric Power: combinations of generators to meet demand
%  - It is possible to skip earlier optional parameters using [] such as:
%       CombinWithLimits(space_per_item, total_space, [], [], item_min)
%    to provide an item_min without specifying an item_max
%
% see also:
%   MultiInvState
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-15 20:47  BryanP      original concept code
%   2  2010-05-16 01:45  BryanP      completed versions 1:3 and tested
%   3  2010-05-16 21:30  BryanP      Added #4: n-D matrix (fastest yet)
%   4  2010-05-26 07:25  BryanP      renamed ValidCombin, expanded comments
%   5  2010-05-26 10:45  BryanP      Added item max, min. Splitout recursion
%   6  2010-05-26 16:50  BryanP      Added min_space_used, debug item_min
%   7  2010-06-06 01:00  BryanP      Added option to specify rounding
%   8  2011-03-08 17:00  BryanP      Expanded usage notes/help
%   9  2011-03-29 10:10  BryanP      Renamed MultiInvValidCombin to CombinWithLimits
%  10  2011-03-30 00:10  BryanP      Corrected initalization of item_counts
%  11  2011-06-10 12:50  BryanP      added value_list and sum_list outputs
%  12  2011-06-10 15:00  BryanP      Fixed missing states using @ceil
%  13  2012-01-19 07:46  BryanP      Sort small_to_big to ensure all possible states are included (esp with @ceil)
%  14  2012-06-02 21:16  BryanP      Bugfix for errant cases with high values for one item and zeros otherwise.
%                                    ALSO: quit for loop early if nothing more to add 
%  15  2012-06-28 21:36  BryanP      Allow specifying fRoundMin, and defaults to complement fRoundMax 
%  16  2012-07-02 14:50  BryanP      BUGFIX: prevent pre-mature loop break when current set empty, but valid future set  
%  17  2016-07-08 07:30  BryanP      Added ability to return summary of input/settings as "opt" 

    %% ----- Initialize optional parameters -----
    %-- Handle total_space
    % set to unbounded if a different bound has not been passed in,
    if isempty(total_space)
        %use no limit
        total_space = Inf;

        %If we have no information about the inventory size, we have a
        %problem...
        if nargin < 3 || isempty(item_max)
            error('CombinWithLimits:NoLimit','must specify either total_space or item_max')
        end
    end

    if nargin < 6 || isempty(fRoundMax)
        fRoundMax = @floor;
    end

    if nargin < 7 || isempty(fRoundMin)
        switch func2str(fRoundMax)
            case 'ceil'
                fRoundMin = @floor;
            case 'round'
                fRoundMin = @round;
            otherwise
                fRoundMin = @ceil;
        end
    end

    %-- Compute max_inv per product
    % if a stricter bound has not been passed in,
    if nargin < 3 || isempty(item_max)
        %simply use the max that will fit in our space (note vectorized)
        item_max = fRoundMax(total_space./space_per_item);
    else
        %otherwise use the value passed in, or the space constraint,
        %which ever is stricter

        %handle scalar by using for all items
        if length(item_max) == 1;
            item_max = item_max*ones(size(space_per_item));
        else
            %verify that our vector is the right size
            assert(length(item_max) == length(space_per_item), ...
                'item_max must have same dimension as space_per_item, or be scalar');
        end
    end

    %-- Handle min_space_used
    % set to zero if a different bound has not been passed in,
    if nargin < 4 || isempty(min_space_used)
        %use no limit
        min_space_used = 0;

    end

    %-- Compute min_inv per product
    % TODO set to just fill min_space_used if a different bound has not been set
        %If we have no user provided information about the minimum size,
    if nargin < 5 || isempty(item_min)
        item_min = zeros(size(space_per_item));
    else
        %otherwise use the values passed in

        %handle scalar by using the same value for all items
        if length(item_min) == 1;
            item_min = item_min*ones(size(space_per_item));
        else
            %verify that our vector is the right size
            assert(length(item_min) == length(space_per_item), ...
                'item_min must have same dimension as space_per_item, or be scalar');
        end
    end

    %Reorder so we add smallest units first. This ensures we include all
    %combinations that totally fill the space, especially with @ceil
    [~, small_to_big] = sort(space_per_item);
    space_per_item = space_per_item(small_to_big);
    item_max = item_max(small_to_big);
    item_min = item_min(small_to_big);

    %Setup our item count so don't have to recompute for each call
    num_items = length(space_per_item);

    %Call recursion
    item_count_list = CombinWithLimitsRecurse(num_items, space_per_item,...
                        total_space, item_max, min_space_used, item_min, ...
                        fRoundMax, fRoundMin);

    %Convert back to original order
    item_count_list(:, small_to_big) = item_count_list;
    space_per_item(small_to_big) = space_per_item;

    if nargout > 1
        value_list = item_count_list .* (ones(size(item_count_list,1),1) * space_per_item);
        sum_list = sum(value_list, 2);
    end
    
    if nargout > 3
        %Return as used input (and default) options if desired
        opt.space_per_item = space_per_item;
        opt.total_space = otal_space;
        opt.item_max = item_max;
        opt.min_space_used = min_space_used;
        opt.item_min = item_min;
        opt.fRoundMax = fRoundMax;
        opt.fRoundMin = fRoundMin;
    end
end

%% ----- Recursive Helper Function -----
function item_counts = CombinWithLimitsRecurse(num_items, space_per_item, ...
                            space_remain, item_max, min_inv_remain, item_min, ...
                            fRoundMax, fRoundMin)

    %When using @ceil or related functions for fRoundMax, it is actually
    %possible to have space_remain <0, if the caller has just barely
    %exceeded the space limit. In that case, all remaining items should be
    %set to zero. Without this check, we end up with a zero dimensioned
    %vector and the valid previous states are not listed. However this
    %"fix" is only valid if it is OK to have a zero value for num_items.
    if space_remain <= 0
        if all(item_min == 0)
            item_counts = zeros(1, num_items);
        else
            item_counts = zeros(0, num_items);
        end
        return
    end

    %Compute the effective item_max & min for this call (vectorized)
    item_max = min(item_max, fRoundMax(space_remain./space_per_item));

    %Basecase (ends the recursion)
    %Note that we know that min & max must be scalars in this case
    %
    % Here is where we enforce the minimum space used constraint, too. Note
    % use of either the overall minimum or the pre-defined item_min. When
    % min > max, item_counts will be empty, which elimates associated item
    % counts that are below the total minimum.
    if num_items == 1;
        item_min = max(fRoundMin(min_inv_remain/space_per_item),item_min);
        item_counts = (item_min:item_max)';
        return;
    end

    %initialize the output array
    % Note: use an empty array with the correct number of columns to avoid
    % warnings and ensure forward compatability
    item_counts = zeros(0,num_items);

    %Handle the first item posibilities here, and other posibilities via
    %recursion
    for num_this_item = item_min(1):item_max(1)

        %recursively get sub-combinations
        other_items = CombinWithLimitsRecurse(num_items-1, ...
            space_per_item(2:end), ...
            space_remain - space_per_item(1)*num_this_item, ...
            item_max(2:end), ...
            min_inv_remain - space_per_item(1)*num_this_item, ...
            item_min(2:end), ...
            fRoundMax, fRoundMin);
        new_rows = size(other_items,1);

        %Note: can't break out early b/c may have impossible combinations
        %now, yet still be possible later.
        if new_rows == 0 && isempty(item_counts)
            continue
        end

        %now add these new combinations to our list of item counts
        %Note: that we prepend the count of the cur_item and num_this item
        %to the new other_items
        item_counts = [ item_counts
                        [ones(new_rows,1)*num_this_item  other_items]
                      ]; %#ok<AGROW>
    end
end %CombinWithLimitsRecurse Function
