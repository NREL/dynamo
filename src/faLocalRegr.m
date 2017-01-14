classdef faLocalRegr < faLocalAvg
%faLOCALREGR local regression function approximation for approx. dynamic programming
%
% Stores a scattered set of n-dimensional known point-value pairs and
% provides function value approximation at other arbitrary points based on
% the linear regression of a neighborhood of nearby points.
%
%USAGE:
% Intitialization:
%   my_approx = faLocalRegr(points, values, options...);
% Adding additional Points:
%   my_approx = update(my_approx, added_points, added_values);
% Approximating:
%   approx_values = approx(my_approx, approx_points, options...);
% Plotting current points and approximation
%   plot(my_approx)
%
%OPTIONS:
% The neighborhood can be specified (ApproxNeighbor) as either all points
% with a given normalized radius (<1) or a specified number of points
% (>=1). This approximation can also be weighted by distance
% (UseDistWeight). With radius neighborhoods, the neighborhood can be
% automatically expanded (AutoExpand) to be a minium specified number
% (PtsOverDim) of points greater than the number of point dimensions.
%
% Points may also be optionally merged together when they lie within a
% specified radius (specify MergeRadius at construction or as a
% parameter to the update function)
%
% The options described above are specifid as string value parameter pairs
% during construction or approximation such as:
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
% before use. In the ADP toolbox these are contined in the kdtreeVER#
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
% See also faLocalAvg, faInterp, faThinPlate
%
% originally by Bryan Palmintier 2012

%% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2012-03-26 09:25  BryanP      Adapted from faThinPlat v2
%   2  2012-03-26 14:50  BryanP      Ball or k-pt ApproxNeighbors & Distance weighting
%   3  2012-03-28 03:00  BryanP      Store ApproxNeighbor Value and DistWeight to match last approx
%   4  2012-03-28 12:30  BryanP      Updated class documentation
%   5  2012-04           BryanP      Added point smoothing
%   6  2012-04-07 13:50  BryanP      Avoid most rank deficient & poorly conditioned problems
%   7  2012-04-21 16:50  BryanP      (re)added autoexpand for ball neighborhoods
%   8  2012-05-06 23:25  BryanP      Extracted faLoaclAvg as super class
%   9  2012-05-08 12:05  BryanP      BUGFIX: fix no points stored and options not used by constructor
%  10  2012-05-14 12:30  BryanP      Rename SmoothRadius to MergeRadius for clarification
%  11  2012-06-20 00:40  BryanP      NEW: MaxRadius, UPDATE: distance weighting
%  12  2012-06-20 01:20  BryanP      Pass approx options to do_approx through obj
%  13  2012-06-20        BryanP      BUGFIX: return vector of values from approx with only one stored sample
%  14  2016-11-10 13:05  BryanP      Expose sampling config for user to edit 


    %Additional properties
    properties
        AutoExpand;         % With neighbor < 1, expand to # of points if needed
        PtsOverDim;         % # points > dimensions for auto-expand (min 1)
    end

    %Internal properties

    methods
        %% ======= Constructor =====
        % Largely uses faLocalAvg, except adding a few extra property defaults
        function obj = faLocalRegr(pts, vals, varargin)

            defaults = {
                        'AutoExpand'        false   % With neighbor < 1, expand to # of points if needed
                        'PtsOverDim'        2       % # points > dimensions for auto-expand (min 1)
                       };
            % Allow for empty contstructor call
            if nargin < 1
                pts = [];
            end
            if nargin < 2
                vals = [];
            end
            if nargin < 3
                varargin = {};
            end

            %Call super class constructor (which sets up the associated
            %defaults)
            obj = obj@faLocalAvg(pts, vals, varargin{:});

            %--Add our approximation specific defaults
            opt = DefaultFields(varargin, defaults);

            %Only set properties that are found in defaults (to prevent
            %altering core functionality)
            for o = 1:size(defaults,1)
                obj.(defaults{o,1}) = opt.(defaults{o,1});
            end

        end

        %% ======= Standard FuncApprox functions =====

        %Note: rely on faLocalAvg for:
        %   update
        %   approx
        %   raw

        function plot(obj, varargin)
            %Disable auto expand for plotting
            old_auto_exp = obj.AutoExpand;
            old_pts_over_d = obj.PtsOverDim;
            obj.AutoExpand = false;
            obj.PtsOverDim = 1;

            plot@faLocalAvg(obj, varargin{:})
            obj.AutoExpand = old_auto_exp;
            obj.PtsOverDim = old_pts_over_d;
        end


    end

    %% HIDDEN METHODS
    methods (Access = protected)
        %% ===== Approximation Specific Functions

        %Note: rely on faLocalAvg for:
        %   build_func

        % ----- do_approx
        % Perform the local regression
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

            % setup kdtree query function as either ball defined by
            % (normalized) radius or a number of nearest points
            if obj.ApproxNeighbor < 1
                idx_lookup_fn = @kdtree_ball_query;
                % for a good regression, we need at least one point point per
                % dimension plus the intercep
                min_num_pts = obj.N_PtDim + obj.PtsOverDim;
            else
                idx_lookup_fn = @kdtree_k_nearest_neighbors;
                if obj.UseDistWeight && obj.ApproxNeighbor < obj.N_PtDim + 1
                    % must have enough points when merging
                    warning('ADP:FuncApprox:LocalRegr:ExpandNeighboor', ...
                            'Distance weighting requires at least N_PtDim + 1 points. Expanding neighborhood from %d to %d', ...
                            obj.ApproxNeighbor, obj.N_PtDim + 1)
                    obj.ApproxNeighbor = obj.N_PtDim + 1;
                else
                    % otherwise allow user to choose fewer points for
                    % nearest neighbor, etc.
                end
                min_num_pts = obj.ApproxNeighbor;
            end

            %Do approximation for each point
            % Note: parfor runs slower b/c of parallel overhead
            for p = 1:size(out_vals,1)
                % find neighborhood of points
                % Note: no need to unnormalize, b/c we have the
                % original, non-normalized point list
                [near_idx, near_dist] = idx_lookup_fn(obj.Func, norm_out_pts(p,:), obj.ApproxNeighbor);

                % remove points that are outside of the Max Radius
                if obj.ApproxNeighbor >= 1 && not(isempty(obj.MaxRadius))
                    too_far_idx = (near_dist > obj.MaxRadius);
                    near_idx(too_far_idx) = [];
                    near_dist(too_far_idx) = [];
                end

                %check that we have enough points for the regression
                if length(near_idx) <  min_num_pts
                    if obj.ApproxNeighbor < 1 && obj.AutoExpand
                        [near_idx, near_dist] = ...
                            kdtree_k_nearest_neighbors(obj.Func, norm_out_pts(p,:), min_num_pts);
                    else
                        %if not, return NaN
                        out_vals(p) = NaN;
                        continue
                    end
                end

                %-- Linear multi-regression
                % extract the neighborhood values
                near_vals = obj.StoreVals(near_idx, :);

                % and the corresponding neighboorhood points
                % Note: tack on a column of ones to compute the intercept
                regr_matrix = horzcat( ones(length(near_idx),1), ...
                                    obj.StorePts(near_idx, :));

                % Remove (obviously) rank deficient columns
                % Note: the MATLAB function rank is more accurate but
                % slow (b/c it uses svd). Here we only catch blantently
                % rank deficient columns with identical values
                valid_cols = (min(regr_matrix(:,2:end)) ~= max(regr_matrix(:,2:end)));

                if any(valid_cols)
                    %include the first column, too
                    valid_cols = [true valid_cols]; %#ok<AGROW> b/c we replace valid_cols above

                    %Check for singular matrices (only an issue if square)
                    if nnz(valid_cols) == length(near_vals)
                        singular_test = rcond(regr_matrix(:,valid_cols));
                        if singular_test < 1e-8
                            %Try removing intercept term
                            valid_cols(1)=false;
                            warning('ADP:FuncApprox:NearSingular', ...
                                'Nearly singular matrix, trying without intercept (RCOND=%g)', ...
                                singular_test);
                        end
                    end

                    % perform the regression
                    if obj.UseDistWeight
                        % Note: distWeight inherited from faLocalAvg
                        near_weights = obj.distWeight(near_dist);
                        %Use lscov for weighted regression
                        local_regr = lscov(regr_matrix(:,valid_cols), near_vals, near_weights);
                    else
                        %Or simply use the backslash operator for basic
                        %regression
                        local_regr = regr_matrix(:,valid_cols)\near_vals;
                    end
                    % extract functional form
                    % if we included the intercept term
                    if valid_cols(1)
                        intercept = local_regr(1);
                        gradient = local_regr(2:end);
                    else
                        intercept = 0;
                        gradient = local_regr;
                    end

                    % now find the value for our point
                    out_vals(p) = intercept + out_pts(p,valid_cols(2:end)) * gradient;
                else
                    %if no valid columns, simply take an average
                    out_vals(p) = mean(near_vals);
                    if size(near_vals,1) > 1
                        warning('ADP:FuncApprox:IdenticalPts', ...
                            '%d identical points found. Using average value', ...
                            size(near_vals,1));
                    end

                end
            end
        end

        %% ===== Helper Methods

        %Note: rely on faLocalAvg for:
        %   normalize
        %   second_nearest_idx

    end

end
