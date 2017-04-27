classdef faLocalAvg < FuncApprox
%faLOCALAVG local average (nearest neighbor) function approximation for approx. dynamic programming
%
% Stores a scattered set of n-dimensional known point-value pairs and
% provides function value approximation at other arbitrary points based on
% the  optionally weighted average of a neighborhood of nearby points.
%
%USAGE:
% Intitialization:
%   my_approx = faLocalRegr(points, values, options...);
% Adding additional Points:
%   my_approx = update(my_approx, added_points, added_values, <MergeRadius>);
% Approximating:
%   approx_values = approx(my_approx, approx_points, <ApproxNeighbor>, <MaxRadius>, <UseDistWeight>);
% Plotting current points and approximation
%   plot(my_approx)
%
%OPTIONS:
% The neighborhood can be specified (ApproxNeighbor) as either all points
% with a given normalized radius (<1) or a specified number of points
% (>=1). This approximation can also be weighted by distance
% (UseDistWeight).
%
% Points may also be optionally merged together when they lie within a
% specified radius (specify MergeRadius at construction or as a
% parameter to the update function)
%
% In the constructor, the options described above are specifid as string
% value parameter pairs during construction or approximation such as:
%  my_approx = faLocalRegr([],[],'ApproxNeighbor', 0.1, 'MergeRadius', 0.001)
% It is also possible to take advantage of MATLAB's cell array expansion
% such as:
%  approx_opt = {'ApproxNeighbor', 0.1, 'MergeRadius', 0.001};
%  my_approx = faLocalRegr([],[], approx_opt{:})
%
%ALGORITHM:
% Uses a kD tree to efficiently find the neighborhood of nearby points.
% Since building the tree is computationally expensive, it is only re-built
% when needed. Such rebuilds are only neccessary the first time approx is
% called and then only if new, unique (not merged) points are added via
% subsequent calls to update. Currently this implementation does not
% support tree-rebalancing, so the entire tree is re-built from scratch
% when needed.
%
%PRIOR TO USING
% Under the hood, the kd_tree implementation uses an open-source kd_tree
% implemented in C++ that must be compiled into platform-specific MEX files
% before use. In the ADP toolbox these are contained in the kdtreeVER#
% subdirectory. Which must also be included in the path. To compile the
% required code use:
%   cd ADP_DIR/kdtree1.2
%   mex kdtree_ball_query.cpp
%   mex kdtree_build.cpp
%   mex kdtree_delete.cpp
%   mex kdtree_k_nearest_neighbors.cpp
%   mex kdtree_nearest_neighbor.cpp
%   (And optionally:) mex kdtree_range_query.cpp
% The above assumes a MATLAB supported compiler has already been setup, use
% doc mex for more information.
%
% See also faLocalRegr, faInterp, faThinPlate
%
% originally by Bryan Palmintier 2012

%% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%  14  2017-04-26 22:55  BryanP      BUGFIX: fix to avoid excess NaN
%  13  2017-04-03 09:18  BryanP      Added autoexpand support
%  12  2016-11-10 13:05  BryanP      Expose sampling config for user to edit 
%  11  2012-06-20        BryanP      BUGFIX: return vector of values from approx with only one stored sample
%  10  2012-06-20 01:20  BryanP      Pass approx options to do_approx through obj
%   9  2012-06-20 00:40  BryanP      NEW: MaxRadius, UPDATE: distance weighting
%   8  2012-06-04 14:30  BryanP      Replace kdtree_nearest_neighbor with kdtree_k_nearest_neighbors(...,1) for improved performance
%   7  2012-05-14 12:30  BryanP      Rename smooth* to merge* for clarification
%   6  2012-05-12 03:35  BryanP      BUGFIX: remove saveobj PLUS: More kdtree handling: copy() and delete()
%   5  2012-05-12 01:45  BryanP      BUGFIX: no points in neighborhood now returns NaN rather than throwing an error
%   4  2012-05-11 23:25  BryanP      Renamed ForceClear/Rebuild as MATLAB special saveobj and loadobj
%   3  2012-05-09 09:35  BryanP      Added ForceClear and ForceRebuild for kD-tree parallel & read/write
%   2  2012-05-07 00:35  BryanP      Added fast special case for neighborhood == 1
%   1  2012-05-06 23:15  BryanP      Adapted from faLocalRegr v7


    properties
        ApproxNeighbor;     %current points to include in approx() and plot(). 0:1 for ball, >1 for # nearest
        MergeRadius;        %(normalized) radius to treat as an identical point and merge
        MaxRadius;          %Largest radius for # neighbor averages
        UseDistWeight;      %Use (normalized) distance when weighting LS fit
        DistFactor;         %distance weighting = 1/(1+(factor*dist)^expon
        DistExpon;          % this is consistant with Martin H, et al 2011
        UseOneAtATime;      %Do not vectorize and cache values (for time comparison)
        AutoExpand;         % With neighbor < 1, expand to # of points if needed
    end

    %Internal properties
    properties (Access=protected)
        %Note: the default FuncApprox property "Func" is used for the kd-tree lookup object
        NPtsMerged;     %Number of points averaged into to the corresponding StorePts
    end

    methods
        %% ======= Constructor =====
        % Largely uses FuncApprox default, but need to change some property defaults
        function obj = faLocalAvg(pts, vals, varargin)

            defaults = {
                        'ApproxNeighbor'    0.1     % <1: normalized radius, >=1: # of points
                        'MergeRadius'       0.001
                        'MaxRadius'         []
                        'UseDistWeight'     false
                        'DistFactor'        20
                        'DistExpon'         1.5
                        'UseOneAtATime'     false
                        'AutoExpand'        true   % With neighbor < 1, expand toat least 1 points if needed
                       };
            % Allow for empty contstructor call
            if nargin < 1
                pts = [];
            end
            if nargin < 2
                vals = [];
            end

            %Call super class constructor
            obj = obj@FuncApprox(pts, vals);

            %--Add our approximation specific defaults
            opt = DefaultFields(varargin, defaults);

            %Only set properties that are found in defaults (to prevent
            %altering core functionality
            for o = 1:size(defaults,1)
                obj.(defaults{o,1}) = opt.(defaults{o,1});
            end

            obj.RefreshIsRequired = true;   %Flag: next approx will require rebuilding approximation
        end

        %% ======= Standard FuncApprox functions =====

        function out_vals = update(obj, new_pts, new_vals, merge_radius)
            if nargin < 4 || isempty(merge_radius)
                %If not specified, set merge radius to stored value
                merge_radius = obj.MergeRadius;
            else
                %If the merging radius has changed, merge-in any cached
                %NewPts using their (old) merging radius
                if (obj.MergeRadius ~= merge_radius) ...
                   && ( (size(obj.NewPts, 1) > 0) ...
                        || (obj.MergeRadius == 0 && obj.RefreshIsRequired && obj.N_StorePts > 0) )
                    obj.build_func();
                end
                obj.MergeRadius = merge_radius;
            end

            %Returning the output values forces an approximation with
            %associated (expensive) tree building, so don't request
            %out_vals unless needed.
            if nargout > 0
                %One at a time is for performance comparison, it forces a
                %tree build (or two if merging required) for each point
                %added
                if obj.UseOneAtATime
                    out_vals = zeros(size(new_vals));
                    for idx = 1:size(new_pts,1)
                        %With outputs to return, we always run an
                        %approximation for each new point
                        out_vals(idx,:) = update@FuncApprox(obj, new_pts(idx,:), new_vals(idx,:), merge_radius);
                    end
                else
                    out_vals = update@FuncApprox(obj, new_pts, new_vals, merge_radius);
                end
            else
                if obj.UseOneAtATime
                    for idx = 1:size(new_pts,1)
                        %With no outputs we force a tre build by calling
                        %build_func after each point is added
                        update@FuncApprox(obj, new_pts(idx,:), new_vals(idx,:), merge_radius);
                        obj.build_func();
                    end
                else
                    update@FuncApprox(obj, new_pts, new_vals, merge_radius);
                end
            end
        end

        function out_vals = approx(obj, out_pts, approx_neighbor, max_radius, dist_weight)
            if nargin >= 3 && not(isempty(approx_neighbor))
                %Update neighborhood if specified
                obj.ApproxNeighbor = approx_neighbor;
            end
            if nargin >= 4
                % Note: empty value is OK for max_radius, it implies
                % ignoring the dist weight
                %
                %Update max radius if specified
                obj.MaxRadius = max_radius;
            end
            if nargin >= 5 && not(isempty(dist_weight))
                %Update dist_weight if specified
                obj.UseDistWeight = dist_weight;
            end

            %After setting our defaults, call our superclass approximation
            %framework. This will in turn call do_approx for the heavy
            %lifting. Options are now passed through the object
            out_vals = approx@FuncApprox(obj, out_pts);
        end

        %Note: rely on default from FuncApprox for:
        %   plot

        function raw_data = raw(obj)
            % Call constructor to include most properites including public
            % and standard hidden properties
            raw_data = raw@FuncApprox(obj);

            % Add on our approximation specific properties
            raw_data.NPtsMerged = obj.NPtsMerged;
        end

    end

    %% HIDDEN METHODS
    methods (Access = protected)
        %% ===== Approximation Specific Functions

        % Construct/update the function approximation
        function build_func(obj)

            % if we only have one point, don't build an approximation
            if obj.N_StorePts == 1
                obj.merge_new_pts();
                obj.NPtsMerged = 1;
                return
            end

            % (1) Setup storage and normalize new points
            % scale to unit hyper-cube (ie [0 1] range for all dimensions
            % Note: normalization includes both old and new points
            norm_store_pts = obj.normalize(obj.StorePts);
            norm_new_pts = obj.normalize(obj.NewPts);
            new_num_pts_merged = ones(size(obj.NewPts, 1), 1);

            % Do merging
            if obj.MergeRadius > 0
                % (2) If there is only one new point, it doesn't need to be
                %merged. But for more new points, attempt to merge any
                %duplicates among them before merging into the stored
                %points

                if size(obj.NewPts, 1) > 1
                    % (3) Build smaller kd-tree of new points
                    try
                        new_pt_tree = kdtree_build(norm_new_pts);
                    catch exception
                        if strcmpi(exception.identifier, 'MATLAB:UndefinedFunction')
                            error('ADP:MissingToolbox', 'kdtree* MEX functions not found, compile them using mex')
                        else
                            rethrow(exception)
                        end
                    end

                    % (4) Merge any new points that are less than the obj.MergeRadius
                    % <a> find nearest among new points
                    nearest_idx = cellfun(@(pt) obj.my_kd_second_nearest(new_pt_tree, pt), ...
                                           num2cell(norm_new_pts,2));
                    % <b> compute distances
                    % Note: realsqrt and realpow assume well formed non-complex
                    nearest_dist = realsqrt(sum(realpow(norm_new_pts - norm_new_pts(nearest_idx,:),2),2));

                    % <c> find items to merge
                    merge_add_map = nearest_dist < obj.MergeRadius;
                    if any(merge_add_map)
                        merge_idx = find(merge_add_map);

                        % <d> merge duplicates
                        to_delete = [];
                        % -i- loop over all remaining duplicates
                        while not(isempty(merge_idx))
                            %-ii- assign the first remaining dup to store results
                            merge_dest = merge_idx(1);

                            %-iii- setup sets of points to merge
                            merge_set = [];
                            add_pts_to_check = merge_dest;

                            %-iv- expand to include all points that are mutually within the merge radius
                            while not(isempty(add_pts_to_check))
                                %a. identify additional points to consider
                                pts_near_add_pts = [];
                                for p = add_pts_to_check;
                                    pts_near_add_pts = vertcat(pts_near_add_pts, ...
                                        kdtree_ball_query(new_pt_tree, norm_new_pts(merge_dest,:), ...
                                            obj.MergeRadius)); %#ok<AGROW>
                                end
                                %b. move all of the points we just checked to the final set
                                merge_set = union(merge_set, add_pts_to_check);
                                %c. see if there are any more point to check
                                add_pts_to_check = setdiff(pts_near_add_pts, merge_set);
                            end

                            %-v- actually do the merging for this point
                            obj.NewPts(merge_dest, :) = mean(obj.NewPts(merge_set, :));
                            norm_new_pts(merge_dest, :) = mean(norm_new_pts(merge_set, :));
                            obj.NewVals(merge_dest, :) = mean(obj.NewVals(merge_set, :));
                            new_num_pts_merged(merge_dest) = sum(new_num_pts_merged(merge_set));

                            %-vi- remove these merged points from the merge_idx list
                            merge_idx = setdiff(merge_idx, merge_set);

                            %-vii- and flag the points (except the merged result) for deletion
                            merge_set = setdiff(merge_set, merge_dest);
                            to_delete = vertcat(to_delete, merge_set); %#ok<AGROW>

                        end
                        % <e> remove duplicates from the new point set
                        obj.NewPts(to_delete, :) = [];
                        norm_new_pts(to_delete, :) = [];
                        obj.NewVals(to_delete, :) = [];
                        new_num_pts_merged(to_delete) = [];

                    % <f> Check for the special case were there is no
                    % existing obj.Func kd-tree AND there were no points to
                    % merge. In this situation, we can use our new tree as
                    % the main tree, and copy over the pts, vals, and merge
                    % count
                    elseif isempty(obj.Func)
                        obj.Func = new_pt_tree;
                        obj.merge_new_pts();
                        obj.NPtsMerged = new_num_pts_merged;
                        return
                    end
                    % <g> clean up new point kd tree.
                    kdtree_delete(new_pt_tree);
                    clear new_pt_tree   %clear this to prevent accidentally using the bad pointer
                end


                % (5) Now check the distances between the (reduced) new set
                % of points and the existing tree (if it exists)
                if not(isempty(obj.Func)) && not(isempty(norm_new_pts))
                    %Note: surprisingly, looping over k_nearest is MUCH,
                    %MUCH faster than calling nearest. In a test with 10000
                    %points, k_nearest was ~50x faster!
                    n_new_pts = size(norm_new_pts, 1);
                    nearest_idx = zeros(n_new_pts,1);
                    for pt = 1:n_new_pts
                        nearest_idx(pt) = kdtree_k_nearest_neighbors(obj.Func, norm_new_pts(pt,:),1);
                    end

                    % (6) Merge any that are less than the merge_radius
                    % <a> compute distances
                    % Note: realsqrt and realpow assume well formed non-complex
                    nearest_dist = realsqrt(sum(realpow(norm_new_pts - norm_store_pts(nearest_idx,:),2),2));

                    % <b> find items to merge
                    merge_add_map = nearest_dist < obj.MergeRadius;
                    if any(merge_add_map)
                        merge_store_idx = nearest_idx(merge_add_map);

                        % <c> average in these points, weighted based on #new vs old
                        %-i- first compute the resulting total merge count
                        obj.NPtsMerged(merge_store_idx) = obj.NPtsMerged(merge_store_idx) + new_num_pts_merged(merge_add_map);
                        %-ii- use this to assign weights
                        m_weight = new_num_pts_merged(merge_add_map)...
                            ./ obj.NPtsMerged(merge_store_idx);
                        %-iii- and merge the values and points
                        obj.StorePts(merge_store_idx) = ...
                             (1-m_weight) .* obj.StorePts(merge_store_idx) ...
                             + m_weight .* obj.NewPts(merge_add_map);
                        obj.StoreVals(merge_store_idx) = ...
                             (1-m_weight) .* obj.StoreVals(merge_store_idx) ...
                             + m_weight .* obj.NewVals(merge_add_map);
                        %-iv- delete the merged points from the new list
                        obj.NewPts(merge_add_map, :) = [];
                        norm_new_pts(merge_add_map, :) = [];
                        obj.NewVals(merge_add_map, :) = [];
                        new_num_pts_merged(merge_add_map) = [];
                    end
                end
            end


            %-- build the lookup tree
            %free memory (at MEX/C level) since we will no longer need our
            %old look up table.
            if not(isempty(obj.Func))
                kdtree_delete(obj.Func);
            end

            %Merge old (store) and new point lists
            obj.merge_new_pts();
            norm_store_pts = vertcat(norm_store_pts, norm_new_pts);
            obj.NPtsMerged = vertcat(obj.NPtsMerged, new_num_pts_merged);

            % construct the kd tree
            try
                obj.Func = kdtree_build(norm_store_pts);
            catch exception
                if strcmpi(exception.identifier, 'MATLAB:UndefinedFunction')
                        error('ADP:MissingToolbox', 'kdtree* MEX functions not found, compile them using mex')
                else
                    rethrow(exception)
                end
            end
        end

        % -------------
        %   do_approx
        % -------------
        % Perform the lookup and local average
        function out_vals = do_approx(obj, out_pts)

            n_out_pts = size(out_pts, 1);

            % if nothing stored, return a vector of NaNs
            if obj.N_StorePts == 0
                out_vals = NaN(n_out_pts,1);
                return
            % if we only have one point, that is our best and only guess
            % and the idea of a normalized neighborhood is also
            % meaningless, so simply return our one value
            elseif obj.N_StorePts == 1
                out_vals = repmat(obj.StoreVals,n_out_pts,1);
                return
            end

            %Initialize output storage
            out_vals = zeros(n_out_pts, 1);

            % normalize our point set to unit hypercube [0 1]
            norm_out_pts = obj.normalize(out_pts);

            % handle special case of nearest neighbor
            if obj.ApproxNeighbor == 1
                near_idx = zeros(n_out_pts,1);
                near_dist = zeros(n_out_pts,1);
                for pt = 1:n_out_pts
                    %Note k_nearest_neighbors is much faster (4-10+x) than
                    %nearest_neighbor. Not sure why.
                    [near_idx(pt), near_dist(pt)] = kdtree_k_nearest_neighbors(obj.Func, norm_out_pts(pt,:),1);
                end

                %Retreive associated values
                out_vals = obj.StoreVals(near_idx, :);
                %Flag those that are too far away
                if not(isempty(obj.MaxRadius))
                    too_far_idx = (near_dist > obj.MaxRadius);
                    out_vals(too_far_idx) = NaN;
                end
            else

                % setup kdtree query function as either ball defined by
                % (normalized) radius or a number of nearest points
                if obj.ApproxNeighbor < 1
                    idx_lookup_fn = @kdtree_ball_query;
                else
                    idx_lookup_fn = @kdtree_k_nearest_neighbors;
                end

                %Do approximation for each point
                % Note: parfor runs slower b/c of parallel overhead
                for p = 1:size(out_vals,1)
                    % find neighborhood of points
                    % Note: no need to unnormalize, b/c we have the
                    % original, non-normalized point list
                    [near_idx, near_dist] = idx_lookup_fn(obj.Func, norm_out_pts(p,:), obj.ApproxNeighbor);

                    %-- Local average
                    % remove points that are outside of the Max Radius
                    if obj.ApproxNeighbor >= 1 && not(isempty(obj.MaxRadius))
                        too_far_idx = (near_dist > obj.MaxRadius);
                        near_idx(too_far_idx) = [];
                        near_dist(too_far_idx) = [];
                    end

                    %Auto expand to include at least one point if needed
                    if isempty(near_idx)
                        if obj.AutoExpand
                            [near_idx, near_dist] = ...
                                kdtree_k_nearest_neighbors(obj.Func, norm_out_pts(p,:), 1);
                        else
                            %if not, return NaN
                            out_vals(p) = NaN;
                            continue
                        end
                    end
                    
                    % extract the neighborhood values
                    near_vals = obj.StoreVals(near_idx, :);

                    % find the average
                    if obj.UseDistWeight
                        near_weights = obj.distWeight(near_dist);
                        out_vals(p) = sum(near_vals .* near_weights)/sum(near_weights);
                    else
                        %If not distance weighting, simply take average
                        out_vals(p) = mean(near_vals);
                    end
                end
            end
        end

        % -------------
        %   distWeight
        % -------------
        function weight_list = distWeight(obj, dist_list)
        % Compute distance weighting
        %
        % With the default values of DistFactor=2 & DistExpon=2, provides
        % near unity weighting of nearby points that falls off to a factor
        % of 1/5 for points 1 normalized distance away.

            weight_list = 1./(1+(obj.DistFactor*dist_list).^obj.DistExpon);
        end

        %% ===== Helper Methods
        function norm_pts = normalize(obj, pts)
            %(1) find range for each dimension
            %Note: PtRange includes stored and new points
            range = obj.PtRange;

            %(2) compute the scaling factor
            dim_range = diff(range,1,1);

            %(3) identify dimensions with zero extent as non-scalable
            dim_to_scale = (dim_range > 0);

            %(4) subtract off minimum
            % Note: bsxfunction is allows implementing:
            %    norm_pts = pts - repmat(range(1,:),n,1);
            % Without the matrix and time overhead of actually building
            % the large repmat matrix
            norm_pts = bsxfun(@minus, pts, range(1,:));

            %(5) scale appropriate dimensions to [0 1]
            % Note: diff(x,1,1) subtracts the bottom row from the top row
            norm_pts(:, dim_to_scale) = ...
                bsxfun(@rdivide, norm_pts(:, dim_to_scale), dim_range(dim_to_scale));

        end

% Unnormalization is not needed b/c we only use the normalized values in
% building the tree, but use the stored, non-normalized values for lookups
%
%         function out_pts = unnormalize(obj, in_pts)
%             range = obj.PtRange;
%
%             %scale back to full [0 1]
%             % Note: diff(x,1,1) subtracts the bottom row from the top row
%             out_pts = bsxfun(@times, in_pts, diff(range,1,1));
%
%             %add back the minimum
%             % Note: bsxfunction is allows implementing:
%             %    pts = pts + repmat(range(1,:),n,1);
%             % Without the matrix and time overhead of actually building
%             % the large repmat matrix
%             out_pts = bsxfun(@minus, out_pts, range(1,:));
%         end

        function second_nearest_idx = my_kd_second_nearest(~, kd_tree, pt)
        % Find the second nearest point index (listed first). Useful
        % to find distinct nearby points when the search point is in the
        % tree. In this case the closest point is the point itself.
        %
        % ONLY HANDLES ONE POINT AT A TIME
            top_two = kdtree_k_nearest_neighbors(kd_tree, pt, 2);
            second_nearest_idx = top_two(1);
        end

        % Override copyElement method to implement custom copy() from mixin
        % so we can force a tree rebuild for the duplicate object
        function cpObj = copyElement(obj)
            % Make a shallow copy of entire object
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % Rebuild the tree
            cpObj.loadobj(obj);
        end

        % Free our tree memory on handle object delete
        function delete(obj)
            if not(isempty(obj.Func))
                kdtree_delete(obj.Func);
            end
        end

%         function sobj = saveobj(obj)
%         % saveobj Prepare kD-tree to be saved (and for parallel uses)
%         %
%         % IMPORTANT: This will crash MATLAB (and maybe your computer) if
%         % used improperly. It calls a low-level (C) memory freeing routine
%         % and hence will cause a Segmentation Fault/Violation if run from
%         % different memory spaces.
%         %
%         % NOTE: saveobj is automatically called by MATLAB during the save
%         % process. It should normally calling this is NOT required since
%         % memory management is handled automatically by update() and
%         % approx().
%         %
%         % IMPORTANT: Use ONLY when saving an approximation to disk AND not
%         % accessed again from this MATLAB session (or a parallel worker).
%         % Additional access must be preceeded by a call to LoadInit.
%         %
%         % In most parallel computing, this will be NOT neccessary since
%         % the new points from workers will likely be merged into the master
%         % using update(). However if an object is to be passed to a worker
%         % and then later used with the worker's new points directly,
%         % SavePrep must be called twice:
%         %   1) Before parfor code to clear main MATLAB session memory pointers.
%         %   2) At end of parfor code
%         %
%         % Forces deletion of the kD-tree.
%         %
%             sobj = copy@FuncApprox(obj);
%             if not(isempty(sobj.Func))
%                 kdtree_delete(sobj.Func);
%                 sobj.Func = [];
%             end
%          end

    end

    methods (Static)
        function obj = loadobj(obj)
        % loadobj Re-initalize kD-tree after loading from disk (and for parallel use)
        %
        % MUST BE DECLARED STATIC
        %
        % NOTE: loadobj is automatically called by MATLAB during the load
        % process.In normal use, this is NOT required since memory
        % management is handled automatically by update() and approx().
        %
        % NOTE: Use ONLY when loading approximations from disk or for
        % parallel computing. For parallel computing it always must be used:
        %   1) at the start of parfor loop body to intialize the approx on
        %   each worker.
        % In addition if an approximation used by a worker is then later
        % used within the master MATLAB session, LoadInit must be called
        % again:
        %   2) after parfor to re-intialize the local approx (note:
        %      updating new points from workers can happen before or after
        %      rebuilding the tree. But rebuilt MUST happen before approx()
        %
        % IMPORTANT: Will cause a memory leak if used in other situations

            if not(isempty(obj.StorePts))
                % normalize Stored points
                norm_store_pts = obj.normalize(obj.StorePts);

                % construct the kd tree
                try
                    obj.Func = kdtree_build(norm_store_pts);
                catch exception
                    if strcmpi(exception.identifier, 'MATLAB:UndefinedFunction')
                            error('ADP:MissingToolbox', 'kdtree* MEX functions not found, compile them using mex')
                    else
                        rethrow(exception)
                    end
                end
            end
        end


    end

end
