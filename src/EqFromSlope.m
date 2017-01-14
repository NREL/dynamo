function [ intercept, slope ] = EqFromSlope(x, slope, y_ref, x_ref)
%EQFROMSLOPE slope & intercept for piecewise linear functions given slopes
%
%   [ intercept, slope ] = EqFromSlope(x, slope, y_ref, x_ref)
%
% Parameters
%   x       vector of x values for the starts of segments
%   slope   vector of slopes for each segment (same size as x)
%   y_ref   (optional) reference y value to fix function vertical location,
%            default = 0. For x_ref = 0, this is the y intercept
%   x_ref   (optional) x value at which to match y_ref, default = x_min
%
% Returns
%   slope   vector of slopes
%   incpt   vector of y intercepts
%
% See also: FunFromSlope, InterpFromSlope, MaxFromSlope

% HISTORY
% ver     date    time        who      changes made
% ---  ---------- -----  ------------- ---------------------------------------
%   1  2011-09-25 19:30  BryanP        Adapted from FunFromSlope v6

%Handle defaults
if nargin < 3
    y_ref = [];
end

if nargin < 4
    x_ref = [];
end

y = FunFromSlope(x, slope, [], y_ref, x_ref);

intercept = y - slope .* x;
end
