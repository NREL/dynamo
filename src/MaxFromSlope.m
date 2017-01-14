function [ x, y_max ] = MaxFromSlope(x_in, slope, x_max, y_ref, x_ref)
%MAXFROMSLOPE finds the maximum for concave piecewise linear functions given slopes
%
%   [ x, y_max ] = MaxFromSlope(x_in, slope, x_max, y_ref, x_ref)
%
% Parameters
%   x_in       vector of x values for the starts of segments
%   slope   vector of slopes for each segment (same size as x)
%   x_max   (optional) parameter to fix the max x value to return a
%            function value for. Extends or truncates as needed
%   y_ref   (optional) reference y value to fix function vertical location,
%            default = 0. For x_ref = 0, this is the y intercept
%   x_ref   (optional) x value at which to match y_ref, default = x_min
%
% Returns
%   x       the x position of the function's max (argmax)
%   y_max   corresponding maximum function value
%
% IMPORTANT: the approximation assumes concavity (slopes monotonically
% decreasing), for non-concave functions, the first local maximum is
% returned in most cases, unless there are two discontinous regions of zero
% slope in which case the results are undefined
%
% Notes:
%  - x_max is optional, but restricts the range of x values returned to be
%    finite. Otherwise if no local max is found with the range of x_in, Inf
%    is returned for x.
%  - For zero slope segments (which correspond to the max) the mid-point x
%    value is returned.
%
% See also: InterpFromSlope, FunFromSlope

% HISTORY
% ver     date    time        who      changes made
% ---  ---------- -----  ------------- ---------------------------------------
%   1  2010-11-04 20:30  BryanP        Adapted from FunFromSlope v6

%=== Compute X value of max ===
%Handle case where x_max is specified
if nargin > 2 && not(isempty(x_max))
    % by truncating the range of slopes to match the range specified by
    % x_max
    slope = slope(x_in<x_max);
else
    %Otherwise assume an infinte x_max
    x_max = Inf;
end

%Now explicitly add x_max
x_in = [x_in(x_in<x_max), x_max];


% Now find the x value of the max taking advantage of the fact that the
% slopes are the first derivative and hence a zero slope (or change in
% sign) marks the max
k_x = find(slope <= 0, 1, 'first');
x = x_in(k_x);

% Handle the case where negative slopes are never found (ie function
% continues increasing toward infinity)
if isempty(x)
    x = x_max;
end


%Handle the special case of a zero slope, by returning the mid-point of the
%zero slope region
if slope(k_x) == 0;
    k_end_flat = find(slope == 0, 1, 'last');
    %Note: since we explicitly tacked on x_max to the end of x_in, but
    %truncated slopes, we know k_end_flat+1 will be a valid index
    x = mean([x, x_in(k_end_flat+1)]);
end


%if just need the x value, skip the y value calculations
if nargout == 1
    return
end

%=== Compute Y value of max ===
%Handle defaults
if nargin < 4 || isempty(y_ref)
    y_ref = 0;
end

if nargin < 5 || isempty(x_ref)
    x_ref = min(x_in);
end

y_max = InterpFromSlope(x_in, slope, x, y_ref, x_ref);

end
