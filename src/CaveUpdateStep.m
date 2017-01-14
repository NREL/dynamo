function [ x, slopes, y, smooth_interval, smooth_k_range, sample_xy] = ...
        CaveUpdateStep(s, s_slopes, u, v, limits, alpha, smooth)
%CAVEUPDATESTEP Single update of Convex Adaptive Value Estimation (ESD.862)
%
%   [ x, slopes, y, smooth_interval, smooth_k_range, sample_xy ] =
%        CaveUpdateStep(s, s_slopes, x, slopes, limits, alpha, smooth)
%
% INPUTS, all required, (Godfrey & Powell name in parens):
%   s          x value of sample (s)
%   s_slopes   1x2 vector of sample slopes below & above s ([pi-, pi+])
%   x          existing x values in estimate (u)
%   slopes     corresponding slopes in estimate (v)
%   limits     [minimum maximum] range for x... this just limits the range
%               of values to approximate over. It is useful to restrict
%               things to avoid undefined values such as log of zero.
%   alpha      smoothing parameter
%   smooth     1x2 vector with minimum smoothing range above&below ([?-,?+])
%
% OUTPUTS
%   x          vector of new x values, these correspond to the start of the interval
%   slopes     vector of coresponding slopes
% ADDITIONAL OUTPUTS
%	y          corresponding function y value (Note: this is an extra
%	           computation using FunFromSlope so not a time saver).
%   smooth_interval
%              actual smoothing interval (over x)
%   smooth_k_range
%              coresponding indexes within the x vector
%   sample_xy  2x1 (X, Y) pair for sample
%
% Special case handling:
%  - out of range (not within limits) are ignored. x_out = x_in and
%    slopes_out = slopes in. Other returns are empty arrays
%  - when the sample x = xmin, only the upper (piplus) slope is used to
%    avoid meaningless non-convexities
%  - non-convex samples slopes are averaged and a warning is issued if the
%    non-convexity is larger than round-off error; however if the first
%    sample (when slopes = [0]) is non-convex (higher slope, piplus, > 0),
%    the higher slope (piplus) is used over the entire range and no error
%    is issued.
%
% Reference:
%  Godfrey, G. A., & Powell, W. B. (2001). An Adaptive, Distribution-Free
%    Algorithm for the Newsvendor Problem with Censored Demands, with
%    Applications to Inventory and Distribution. Management Science, 47(8),
%    1101-1112.
%
% see also:
%   CaveUpdate, FunFromSlope
%
% originally by Bryan Palmintier 2010
%   rewritten from initial attempts by Sarvee D & Mort Webster

% TODO:
%  - rename alpha to step
%  - streamline when using zero smoothing interval
%  - rework alpha handling to work for a local area rather than across the
%    entire function
%  - Convert to an object oriented paradigm
%  - provide defaults for many options
%
% Implementation Note:
%  This version (3+) uses a modified version of the algorithm that first
%  smooths the minimum smoothing interval and then extends it if needed.
%
% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-04-18 01:14  BryanP      Split from CaveUpdate ver 6
%   2  2010-04-18 02:00  BryanP      Ooops... added the smoothing back
%   3  2010-04-20 11:00  BryanP      Converted to Smooth first then extend
%   4  2010-04-21 07:50  BryanP      Updated to:
%                                     - remove redundant points via unique
%                                     - optionally compute the function values
%   5  2010-08-04 15:00  BryanP      Added warning for non-convex sample slopes
%   6  2010-08-10 22:21  BryanP       a) Added fix for sample = x_min
%                                     b) Fixed upper smoothing interval epsplus
%                                     c) replaced unused "junk" returns with ~
%                                     d) removed unused u_ub and u_lb calcs
%   7  2010-08-11 11:30  BryanP      Added convexity tolerance to avoid excess
%                                    non-convex warnings caused by roundoff
%   8  2010-08-11 22:15  BryanP      Reverted to "junk" for compatability
%                                    with older MATLAB versions.
%   9  2010-11-02 23:55  BryanP      Added special case to use piplus for
%                                    entire range on first pass
%  10  2010-11-02 23:55  BryanP      Refined first pass handling to only
%                                    extend piplus to zero for non-convex
%                                    samples to avoid averaging against
%                                    in-correct negative slopes.
%  11  2011-03-10 12:28  BryanP      NaN check & correct for slopes

%Uncommmonly changed options
tol = 1e-10;    %Tolerance for difference between slopes piminus & piplus

%--- Parse input parameters
% Partially adapted to Godfrey & Powell nomenclature
piminus = s_slopes(1);
piplus = s_slopes(2);

epsminus = smooth(1);
epsplus = smooth(2);

x_min = limits(1);
x_max = limits(2);

%set flag if we need to compute the function's y values or sample xy...
compute_y =  nargout > 2;

%Handle out of range samples by ignoring them
if (s > x_max) || (s < x_min)
    warning('CAVEStep:SampleOutOfRange', ...
        'Ignoring out of range sample (sample_x=%g not in limits = [%g %g]',s, x_min, x_max)
    x = u;
    slopes = v;
    if compute_y
        smooth_interval = [];
        smooth_k_range =  [];
        sample_xy = [];
    end
    return
end

%Handle NaNs by attempting to use the other slope
if any(isnan(s_slopes))
	warning('CAVEStep:NaNSlope', ...
	'NaN Sample slopes encountered, attempting to use remaining valid slope [%g %g]',...
        piminus, piplus)
    if isnan(piminus)
        if isnan(piplus)
            error('CAVEStep:TwoNaNSlopes', ...
                'Failed. Both slopes are NaN [%g %g]', piminus, piplus)
        else
            %piplus is a valid number so use it for both
            piminus=piplus;
        end
    else
        %piminus is a valid number so use it for both
        piplus=piminus;
    end
end

%Handle values at xmin by using the slope above (piplus) for both values.
%This avoids non-convexities when the slope below (piminus) is undefined or
%set to zero
if s == x_min
    piminus = piplus;
end

%Handle non-convex sample slopes by taking average (prevents odd
%re-ordering of x values, and ensures a truly convex result
if piminus < piplus
    % Handle special case of the first non-zero sample with a zero lower
    % slope which is non-convex by definition if piplus is positive. This
    % situation often arises when starting with a flat value function
    % approximation. In this case a better convex estimate is to simply take
    % piplus over the entire region, rather than using an average. This is
    % particularly true with rapidly decreasing stepsizes where the
    % combination decreased stepsizes and a poor (averaged, often with
    % zero) initial guess may take many samples to correct the slopes
    % corresponding to higher x values. End result i better initial
    % estimates early in approximations.
    %
    % Note: we only do this for non-convex cases (typically with a zero
    % lower slope, piminus. If the slope pair is convex (typically implying
    % a negative piplus) we don't want to propigate the negative value back
    % to zero, because rapidly decresing step sizes will have a difficult
    % time correcting the situation, resulting in persistantly incorrect
    % approximations in practice.
    if length(v) == 1 && v == 0;
        x = u;
        slopes = piplus;
        if compute_y
            smooth_interval = limits;
            smooth_k_range =  [1 1];
            sample_xy = [s, s*slopes];
        end
        return
    end

    %Otherwise, use the average
    piavg = mean([piminus piplus]);

    %provide warning only if non-convexity exceeds tolerance (suppresses
    %superflous warnings for round-off errors
%    if piminus + tol < piplus
%         warning('CAVEStep:NonConvexSampleSlopes', ...
%             'Sample slope below larger than above (not convex): [%g < %g], using average=%g',...
%              piminus, piplus, piavg)
%    end
    piminus = piavg;
    piplus = piavg;
end

% ---- Insert Sample
%-- insert our sample if needed
% Note: this makes it much easier to limit the search ranges if we need
% to extend the smoothing interval.
k_s = find(u<=s,1,'last');
%should always be found, b/c we verified that s is in range
if u(k_s) ~= s
    %if it doesn't exist, insert it, and use old slope for segment
    %Note: k_s will match the k less than the sample so need to add 1
    k_s=k_s+1;
    [u, v] = CaveInsert(k_s, s, u, v);
end

%% ---- Identify Minimum Smoothing Interval
% -- Lower Bound (of minimum smoothing interval)
% limit the lower bound to within the limits
u_lb = max(x_min,s-epsminus);
% insert the lower bound breakpoint if it doesn't already exist by
%  first looking for the highest current value below it
k_lb = find(u<=u_lb,1,'last');
if isempty(k_lb)
    %We should only get here if the user passed us a set of u's that
    %doesn't extend to x_min. This should never happen
    %but if it does, give a warning...
    warning('CAVE:MissingXmin', ...
        'Passed x does not extend to lower limit (it should), adding')
    % and add the lower bound for the user
    u = [x_min, u];
    v = [v(1), v];
    k_lb = 1;
end

%if it doesn't exist, insert it, and use old slope
if u(k_lb) ~= u_lb
    %Note: k_lb will match the k less than the sample so need to add 1
    k_lb=k_lb+1;
    [u, v] = CaveInsert(k_lb, u_lb, u, v);

    %don't forget that we (probably) inserted something before our sample
    %unless the lower bound and sample are in the same place
    if u_lb ~= s
        %so need to update our pointer index
        k_s=k_s+1;
    end
end

% -- Upper Bound (of minimum smoothing interval)
% limit the upper bound to within the limits
u_ub = min(x_max,s+epsplus);
% insert the upper bound breakpoint if it doesn't already exist by
%  first looking for the lowest current value above it
k_ub = find(u>=u_ub,1,'first');

%if we can't find a u>u_ub, then we should smooth to the end of our
%interval
if isempty(k_ub)
    %We should only get here if the upper bound is higher than the
    %currently highest u value, which is fine,
%     %so simply add the upper bound for the user
%     u(end+1) = u_ub;
%     % and take the last current slope
%     v(end+1) = v(end);
%
%     %and be sure to update our indices
    %It's tempting to add this upper bound to our point set in the
    %approximations, but this will actually be redundant, since the slope
    %will be the same as that from the last current value. Therefor, we
    %simply put the upper bound at our highest u, so that k_ub is a valid
    %index.
    k_ub=length(u);

    %choosing to update u_ub or not is a stylistic choice... here we opt
    %not to.
%    u_ub=u(k_ub);

%otherwise, if it doesn't exist, insert it, and use old slope
elseif u(k_ub) ~= u_ub
    [u, v] = CaveInsert(k_ub, u_ub, u, v);
end

% -- Perform the initial smoothing
%lower half of smoothing interval
v = CaveSmooth(k_lb, k_s-1, piminus, v, alpha);
%upper half of smoothing interval
v = CaveSmooth(k_s, k_ub, piplus, v, alpha);

%% ---- Now check for convexity and perform additional smoothing if needed
% Note: the key to understanding this approach is recognizing that a convex
% line has monotonically decreasing slopes, and so should/must the existing
% vector of slopes passed to this function.

%-- First check below our smoothed range... in this area, we need to
%identify any slopes that are lower than our first smoothed slope, at k_lb
k_minus = find(...
    v(1:k_lb-1) ...
    < CaveSmooth(1, k_lb-1, piminus, v(2:k_lb), alpha), ...
    1,'first');

%if there are any, smooth back to there.
% note: assuming the original slopes were convex, this will ensure the
% final slopes are also convex, because the adjustment factor
% (alpha+piminus...) is the same.
if ~isempty(k_minus)
    v = CaveSmooth(k_minus, k_lb-1, piminus, v, alpha);
    k_lb = k_minus;
    %    u_lb = u(k_minus);    %Only needed for debugging
end

%-- And do the same for above our smoothed range... in this area, we need
%to identify any slopes that are higher than our last smoothed slope,
% at k_ub
k_plus = k_ub + find(...
    CaveSmooth(1, length(v)-k_ub, piplus, v(k_ub:end-1), alpha) ...
    < v(k_ub+1:end), ...
    1,'last');

%if there are any, smooth forward to there.
% note: assuming the original slopes were convex, this will ensure the
% final slopes are also convex, because the adjustment factor
% (alpha+piplus...) is the same.
if ~isempty(k_plus)
    v = CaveSmooth(k_ub+1, k_plus, piplus, v, alpha);
    k_ub = k_plus - 1;
    %    u_ub = u(k_ub);    %Only needed for debugging
end

%% ---- Remove any redundant points
% Explain

% -- First compute the sample xy, if needed, before we possibly remove the
% sample
if compute_y
   y = FunFromSlope(u,v);
   sample_xy = [u(k_s), y(k_s)];
end

% -- Now remove redundant points

% Using MATLAB's unique function
% Notes:
%   - by default, unique sorts in ascending order, we want the reverse
%   - we must specify 'first' to keep the proper u values
[junk, convert2short new_idx] = unique(v, 'first'); %#ok<ASGLU>
convert2short = convert2short(end:-1:1);
new_idx = new_idx(end:-1:1);

%Now update the actual variables, mantaining the correct x,y pairing
v = v(convert2short);
u = u(convert2short);

if compute_y
    y = y(convert2short);

    %and update intervals as needed
    k_lb = new_idx(k_lb);
    k_ub = new_idx(k_ub);
    smooth_interval = [u(k_lb), u(k_ub)];
    smooth_k_range =  [k_lb k_ub];
end

%% -- Outputs
x = u;
slopes = v;
end

function [x, slopes] = CaveInsert(k, x_new, x, slopes)
    %CaveInsert helper function for CaveUpdateStep
    %
    % Inserts new_x into x at k and copies the current segment slope
    %
    % IMPORTANT: Only works for k>=2
    x = [x(1:k-1), x_new, x(k:end)];
    slopes = [slopes(1:k-1), slopes(k-1), slopes(k:end)];
end

function [slopes] = CaveSmooth(k_min, k_max, new_slope, slopes, alpha)
    %CaveSmooth helper function for CaveUpdateStep
    %
    % Smooths the defined interval
    slopes(k_min:k_max) = alpha*new_slope + (1-alpha)*slopes(k_min:k_max);
end
