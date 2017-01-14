function [ y, x_out ] = FunFromSlope(x, slope, x_max, y_ref, x_ref)
%FUNFROMSLOPE computes y values for piecewise linear functions given slopes
%
%   [ y, x_out ] = FunFromSlope(x, slope, x_max, y_ref, x_ref)
%
% Parameters
%   x       vector of x values for the starts of segments
%   slope   vector of slopes for each segment (same size as x)
%   x_max   (optional) parameter to fix the max x value to return a
%            function value for. Extends or truncates as needed
%   y_ref   (optional) reference y value to fix function vertical location,
%            default = 0. For x_ref = 0, this is the y intercept
%   x_ref   (optional) x value at which to match y_ref, default = x_min
%
% Returns
%   y       vector of function y values
%   x_out   vector of x values corresponding to y
%
% Note that x_max is optional, but allows the user to force the function to
% extend beyond the current range of x.
%
% See also: InterpFromSlope, MaxFromSlope, EqFromSlope

% HISTORY
% ver     date    time        who      changes made
% ---  ---------- -----  ------------- ---------------------------------------
%   1  2010-04-17 14:30  BryanP        Initial code
%   2  2010-08-10 16:00  MortW         Added quick fix for intercept
%   3  2010-08-10 22:47  BryanP        Made intercept parameter optional
%   4  2010-08-10 23:59  BryanP        a) added support for arbitrary x&y references
%                                      b) Renamed intercept to y_ref
%   5  2010-08-11 00:53  BryanP        Added support for out of range x_ref
%   6  2010-11-04 16:00  BryanP        Use InterpFromSlope for x_ref calcs

%Handle defaults
if nargin < 4 || isempty(y_ref)
    y_ref = 0;
end

%Handle specified, non-standard x_ref
if nargin >= 5 && not(isempty(x_ref)) && x_ref ~= min(x);
    %In this case we first compute the corresponding function without any
    %offset references and then interpolate to compute the vertical
    %position
    y_baseline = InterpFromSlope(x, slope, x_ref, 0);

    y_ref = y_ref - y_baseline;
end


%Handle x_max, there are two conditions both handled by the code below:
% 1) the range of x extends beyond x_max, in which case x_out and slope
% must be truncated
% 2) the range of x does not reach x_max (or the max value of x after truncation
% is less than x_max), in which case we need to tack on x_max. There is
% actually no need to adjust the slopes in this case, because the final
% slope value is not needed
if nargin >= 3 && not(isempty(x_max))
    x_out = [x(x<x_max) x_max];
    slope = slope(x<x_max);
else
    %If there is no x_max specified, simply use the entire range of x
    x_out = x;
    slope = slope(1:(length(x_out)-1));
end

y = y_ref + [0, cumsum(diff(x_out).*slope)];
end
