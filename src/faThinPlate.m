classdef faThinPlate < FuncApprox
%faTHINPLATE thin plate function approximation for approx. dynamic programming
%
% Uses the MATLAB Curve Fitting Toolbox's thin plate approximation for
% smooth function approximation of scattered 2-D points
%
% See also faLocalRegr, faInterp
%
% originally by Bryan Palmintier 2012

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2012-03-25 15:45  BryanP      Adapted from faDiscrete v7
%   2  2012-03-25 22:15  BryanP      Abstract out core functionality to FuncApprox superclass
%   3  2012-03-26 14:50  BryanP      Hide approximation specific functions


    methods
        %% ======= Constructor =====
        % Largely uses FuncApprox default, but need to change our defaults
        % for some properites
        function obj = faThinPlate(pts, vals, varargin)
            if nargin < 1
                pts = [];
            end
            if nargin < 2
                vals = [];
            end
            %Call super class constructor
            obj = obj@FuncApprox(pts, vals);

            %Update our approximation specific values
            obj.MinPtDim = 2;
            obj.MaxPtDim = 2;
            obj.RefreshIsRequired = true;   %Flag: next approx will require rebuilding approximation
        end

        %% ======= Standard FuncApprox functions =====
        % Use FuncApprox defaults


    end

    %% HIDDEN METHODS
    methods (Access = protected)
        %% ===== Approximation Specific Functions
        function build_func(obj, varargin)
            obj.merge_new_pts();

            try
                %Ergh... the curve fitting toolbox wants points in
                %columns, so we transpose & transpose back later
                obj.Func = tpaps(obj.StorePts', obj.StoreVals');
            catch exception
                if strcmpi(exception.identifier, 'MATLAB:UndefinedFunction')
                    error('ADP:MissingToolbox', 'faThinPlate requires the Curve Fitting Toolbox')
                else
                    rethrow(exception)
                end
            end
        end

        function out_vals = do_approx(obj, out_pts, varargin)
            %Ergh... the curve fitting toolbox wants points in
            %columns, so we transpose & transpose back later
            out_vals = fnval(obj.Func, out_pts')';
        end


    end

end
