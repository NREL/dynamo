function y = InterpFromSlope(x, slope, x_out, y_ref, x_ref)
%INTERPFROMSLOPE specified y values for piecewise linear functions given slopes
%
%   y = InterpFromSlope(x, slope, x_out, y_ref, x_ref)
%
% Parameters
%   x       vector of x values for the starts of segments
%   slope   vector of slopes for each segment (same size as x)
%   x_out   x value(s) for which to compute the corresponding y
%   y_ref   (optional) reference y value to fix function vertical location,
%            default = 0. For x_ref = 0, this is the y intercept
%   x_ref   (optional) x value at which to match y_ref, default = x_min
%
% Returns
%   y       vector of function y values
%
% NOTE: This version of this function is best for a small number of x_out values, because
% it approximates each x_out value independantly. To compute a large number
% of x_out values use FunFromSlope() with interp1q()
%
% TODO: Vectorize for support of larger x_out vectors
%
% See also FunFromSlope, MaxFromSlope

% HISTORY
% ver     date    time        who      changes made
% ---  ---------- -----  ------------- ---------------------------------------
%   1  2010-04-17 14:30  BryanP        Adapted FunFromSlope v5
%   2  2011-01-31 23:40  BryanP        Updated comments

%Handle defaults
if nargin < 4 || isempty(y_ref)
    y_ref = 0;
end

%Handle specified, non-standard x_ref
if nargin >= 5 && not(isempty(x_ref)) && x_ref ~= min(x);
    %In this case we first compute the corresponding function without any
    %offset references and then interpolate to compute the vertical
    %position (uses recursion)
    y_baseline = InterpFromSlope(x, slope, x_ref, 0);

    y_ref = y_ref - y_baseline;

end

%allocate space and match x_out shape
y = zeros(size(x_out));

%approximate each x_out individually
%TODO - Vectorize
for idx = 1: length(x_out)
    this_x = x_out(idx);
    if this_x > x(1)
        y(idx) = y_ref + sum(diff([x(x<this_x) this_x]).*slope(x<this_x));
    elseif this_x == x(1)
        y(idx) = y_ref;
    else
        error('ADPtoolbox:InterpFromSlope:x_out_of_range', ...
            'InterpFromSlope: cannot estimate value for x_out < min(x)')
    end
end
