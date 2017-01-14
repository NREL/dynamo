function [ u, v] = ...
    CaveUpdate(s, model, u, v, limits, alpha_set, smooth_set, sched, visual)
%CAVEUPDATE  Convex Adaptive Value Estimation for multiple samples(ESD.862)
%
%   [ x, slopes] = ...
%    CaveUpdate(samples, @model, x, slopes, limits, alpha, smooth, sched, visual)
%
% inputs, all required, (Godfrey & Powell name in parens):
%   samples    vector of n samples (s)
%   model      function handle to call for estimating slopes at s
%   limits     [minimum maximum] range for x
%   alpha      smoothing factor
%   smooth     corresponding minimum smoothing range ([?-,?+])
%   sched      iterations to start using corresponding alpha and smooth as
%              follows:
%                alpha/smooth  schedule   effect
%                ------------  --------   --------------------------------
%                                 []      Constant alpha(1) & smooth(1)
%                  1/2 x z        1xz     change to new values on sched iteration
%                 scalar/1x2      1xz     reduce values by 50% on sched iterations
%                 scalar/1x2    scalar    reduce by 50% every sched iterations
%   visual     0/omit: no plot, int: show every VISUAL points on graph
%
% Notes:
%  - alpha & smoothing parameter schedule not fully tested
%
% Reference:
%  Godfrey, G. A., & Powell, W. B. (2001). An Adaptive, Distribution-Free
%    Algorithm for the Newsvendor Problem with Censored Demands, with
%    Applications to Inventory and Distribution. Management Science, 47(8),
%    1101-1112.
%
% see also:
%   CaveUpdateStep, CaveTest
%
% originally by Bryan Palmintier 2010
%   rewritten from initial attempts by Sarvee D & Mort Webster

% HISTORY
% ver     date    time        who      changes made
% ---  ---------- -----  ------------- ---------------------------------------
%   0  2009-?            MortW         Initial, not working (?) code
%   1  2010-04-15 AM     SarveeD       Revised code... still not quite
%   2  2010-04-15 PM     SarveeD       Mostly working code
%   3  2010-04-16 13:30  BryanP        Added visual debugging & randperm of obs
%   4  2010-04-17 22:00  BryanP        It works! Completed overhaul of code:
%                                       - vectorized smooth interval calcs
%                                       - handled end cases for k & Q
%                                       - removed hard-coded constants
%   5  2010-04-17 23:00  BryanP        Fixed kplus bug on insert of u_lb
%   6  2010-04-18 01:14  BryanP        Converted to function:
%                                        - core algorithm in CaveUpdateStep
%                                        - problem specific cave_scratch
%   7  2010-04-21 11:00  BryanP        Updated for CaveUpdateStep v4

%% ---- Handle inputs
if nargin < 9
    visual = 0;
end

x_max = limits(2);
smooth_scale = 0.5;

%% Initialize plot if needed
if visual
    [myF,uInit]=FunFromSlope(u,v,x_max);
    %this hold on is for existing figures
    hold on
    plot(uInit,myF)
    %this hold on is for new figures
    hold on
end

%% Initialize smoothing parameters
alpha = alpha_set(1);
smooth = smooth_set(1,:);
%setup next time to adjust the parameters
if isempty(sched)
    smooth_adjust = Inf;
else
    smooth_adjust = sched(1);
end
smooth_idx = 1;

%% Loop through samples
for ctr = 1:length(s)
    %-- Adjust smoothing parameters if needed
    if ctr == smooth_adjust
        smooth_idx = min(smooth_idx+1,max(length(alpha),length(sched)));
        if length(alpha) == 1;
            alpha = alpha*smooth_scale;
            smooth = smooth*smooth_scale;
        else
            alpha = alpha_set(smooth_idx);
            smooth = smooth_set(smooth_idx,:);
        end
        if length(sched) == 1
            smooth_adjust = smooth_adjust + sched;
        else
            smooth_adjust = sched(smooth_idx);
        end

    end

    %-- Get sample slopes
    s_slope = model(s(ctr));

    %-- Run CAVE algorithm & Update Plot if desired
    if mod(ctr,visual) ~= 0
        %if not plotting, save the computation of y values (myF)
        [u, v] = CaveUpdateStep(s(ctr), s_slope, u, v, limits, alpha, smooth);
    else
        %but if plotting, collecte all the data we can
        [u, v, myF, junk_interval, k, sample_xy] = ...
            CaveUpdateStep(s(ctr), s_slope, u, v, limits, alpha, smooth);
        %we have the function... but need to add x_max
        myF = [myF, myF(end)+v(end)*(x_max-u(end))];
        plot_x = [u x_max];

        %highlight smoothing interval
        % Note: add 1 to k_ub to extend to end of interval
        k_lb = k(1);
        k_ub = k(2);
        k_ub = min([k_ub+1,length(myF)]);
        plot(plot_x(k_lb:k_ub),myF(k_lb:k_ub),'Color', 0.8*[1 1 1],'LineWidth', 2)
        %plot function (after highlighting interval, so no written over)
        plot(plot_x,myF)

        %plot sample
        x_s = sample_xy(1);
        y_s = sample_xy(2);
        %highlight sample
        plot(x_s, y_s, 'x')
        text(x_s, y_s, int2str(ctr))
     end
end

if visual
    [myF,plot_x]=FunFromSlope(u,v,x_max);
    plot(plot_x,myF,'r')
    hold off
end
