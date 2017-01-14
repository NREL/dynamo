function [ x, slopes, samples ] = ...
    CaveTest(my_fun, x_range, samples, fun_name, alpha, smooth, sched, ds)
%CaveTest Simplified testing of the Cave Algorithm (ESD.862)
%   Automates the plotting of the sampling, derivatives & plotting of 
%   simple functions using the CAVE algorithm
%
% [ x, slopes, samples ] = ...
%     CaveTest(@my_fun, x_range, samples, fun_name, alpha, smooth, sched, ds)
%
% Inputs: (defualts for optional parameters in parentheses)
%   @my_fun   handle for vector capable value function. CaveTest handles
%             finite differencing
%   x_range   [x_min, x,max]
%   samples   either:
%               - scalar: number of uniform random samples on [x_min+ds,
%                 x_max-ds]
%               - vector of samples to try
%   fun_name  string to use in plot title for function name (func2str(my_fun)
%   alpha     scaling factor... see CaveUpdate() for options (1)
%   smooth    smoothing interval... see CaveUpdate() options ([2 2])
%   sched     smoothing schedule... see CaveUpdate() options (25)
%   ds        the interval to use for finite differences (0.1)
%
% Note: a flat start is used with x=x_min & slope(x_min)=0

% HISTORY
% ver     date    time        who      changes made
% ---  ---------- -----  ------------- ---------------------------------------
%   1  2010-04-20 03:00  BryanP        Adapted from cave_scratch ver 2

%-- Handle input arguments
%optional arguments
if nargin < 4 || isempty(fun_name)
    fun_name = func2str(my_fun);
end

if nargin < 5 || isempty(alpha)
    alpha = 1;
end

if nargin < 6  || isempty(smooth)
    smooth = [2 2];
end

if nargin < 7 || isempty(sched);
    sched = 25;
end

if nargin <8 || isempty(ds);
    ds = 0.1;
end

%split-up inputs as needed
x_min = x_range(1);
x_max = x_range(2);

%handle samples if needed
if length(samples) == 1;
    %Note adjustment by 2ds to ensure we don't sample outside of the range
    samples = rand(1,samples)*(x_max-x_min-2*ds)+x_min+ds;
end

n_samples = length(samples);


%% -- Setup visualizations
%plot actual values
plot_x=linspace(x_min,x_max,100);
plot_y=my_fun(plot_x);

plot(plot_x,plot_y,'g','LineWidth',2)
title(sprintf('%s with flat start & %d random samples',fun_name, n_samples))

%% --Setup and run CAVE
% Initial Approximation
u = x_min;
v = 0;

% Derivative function
my_slopes = @(s) diff(my_fun([s-ds,s,s+ds]))/ds;


%Actually run the algorithm & make pretty pictures
[u,v] = CaveUpdate(samples, my_slopes, u, v, [x_min x_max], ...
                           alpha, smooth, sched, ceil(n_samples/20));


if nargout > 1
    x = u;
    slopes = v;
end
end

