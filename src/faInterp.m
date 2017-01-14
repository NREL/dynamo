classdef faInterp < FuncApprox
%faINTERP Scattered interpolation function approximation for approx. dynamic programming
%
% Uses MATLAB's TriScatteredInterp to interpolate scattered 2-D or 3-D
% points using Delaunay triangulation
%
% Method options: 'natural', 'linear', 'nearest'
%
% See also faLocalRegr, faInterp, faThinPlate
%
% originally by Bryan Palmintier 2012

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2012-03-25 22:30  BryanP      Adapted from faThinPlate v2
%   2  2012-03-26 15:00  BryanP      Allow on-the-fly method changes
%   3  2012-06-20 02:00  BryanP      BUGFIX: use stored method when unspecified
%   4  2012-07-06 17:30  BryanP      Suppress duplicate point avg warning

    %Read only properties
    properties (SetAccess=protected)
        InterpMethod = 'natural';    %Delaunay interpolation method
    end

    methods
        %% ======= Constructor =====
        % Largely uses FuncApprox default, but need to change our defaults
        % for some properties
        function obj = faInterp(pts, vals, method, varargin)
        % obj = faInterp(pts, vals, method)
            if nargin < 1
                pts = [];
            end
            if nargin < 2
                vals = [];
            end
            %Call super class constructor
            obj = obj@FuncApprox(pts, vals);

            if nargin >= 3
                obj.InterpMethod = method;
            end

            %Update our approximation specific values
            obj.MinPtDim = 2;
            obj.MaxPtDim = 3;
            obj.RefreshIsRequired = true;   %Flag: next approx will require rebuilding approximation
        end

        %% ======= Standard FuncApprox functions =====

        % Use FuncApprox defaults
        function out_vals = approx(obj, out_pts, method)
            if nargin < 3 || isempty(method)
                method = obj.InterpMethod;
            else
                obj.InterpMethod = method;
            end

            %if the interpolation method has changed, we need to update our
            %approximation
            if not(isempty(obj.Func)) && not(strcmp(obj.Func.Method, method))
                obj.RefreshIsRequired = true;
            end
            out_vals = approx@FuncApprox(obj, out_pts, method);
        end


    end

    %% HIDDEN METHODS
    methods (Access = protected)
        %% ===== Approximation Specific Functions
        function build_func(obj, method)
            if nargin >=2 && not(isempty(method))
                obj.InterpMethod = method;
            end
            obj.merge_new_pts();

            %Temporarily suppress duplicate point warnings (and save state)
            old_warn_state = warning('off', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId');

            obj.Func = TriScatteredInterp(obj.StorePts, obj.StoreVals, obj.InterpMethod);

            % Restore the state of that duplicate point warning
            warning(old_warn_state)

        end

        function out_vals = do_approx(obj, out_pts, ~)
            out_vals = obj.Func(out_pts);
        end


    end

end
