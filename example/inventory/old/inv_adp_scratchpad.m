% Inventory ADP Scratchpad
%
% Contains example runs of the inventory_adp algorithm for testing and to
% demonstrate usage and functionality.

% HISTORY
% ver     date    time        who      changes made
% ---  ---------- -----  ------------- ---------------------------------------
%   1  2010-11-02 22:30  BryanP        Initial version 
%   2  2010-11-03 02:30  BryanP        Added super simple demonstration call 
%   3  2010-11-03 23:30  BryanP        Tuned medium case & compute error
%   4  2010-11-07 22:30  BryanP        Renamed configurable simple case to
%                                       med and med to big. Added
%   5  2012-02-29 14:00  BryanP        Added summary for "big" iterations

%% ========== A super simple call demonstration ==========
%Note: this uses the same defaults from inventory_adp v10

%IMPORTANT: the value function approximations have changed. They are now
% the true post decision state and hence may not appear to have a finite
% maximum. This is correct according to the ADP forms of Bellmans found in
% Powell (2007). When ordering costs are introduced, the curve is
% effectively tilted such that there is a finite maximum as seen in the
% orders matrix

%options
N = 100;
iter = 1000;

fprintf('\n\nStarting simple (default) call to inventory_ADP\n')
fprintf('  n_periods = %d, iterations = %d\n', N, iter)
tic
[simple_ADP_orders, simple_ADP_results, simple_ADP_values] = inventory_adp(N, iter);
toc

%% ========== The companion traditional DP call ==========
% First run exact dp for comparison
N = 100
disp('Default DP basecase for comparision')
fprintf('  n_periods = %d\n', N)
tic
[simple_DP_orders, simple_DP_values] = inventory3(50, N, 20, 0.05);
toc

%% ========== companion traditional DP call with DP function =======
clear;
clc;

N = 100;

maximum_inventory = 20;

disp('Default DP basecase for comparision using DP function')
fprintf('n_periods = %d\n', N)

parameters.maximum_inventory = maximum_inventory;

poisson_lambda = 20;

parameters.p_demand = poisspdf(0:(maximum_inventory-1), poisson_lambda)';
parameters.p_demand(maximum_inventory+1) = 1 - sum(parameters.p_demand);


parameters.sales_price = 8;
parameters.terminal_unit_value = 0;
parameters.order_cost = -4;
parameters.unit_cost = -2;
parameters.hold_cost = -1;

tic
[DP_policy, DP_cost_to_go] = DP(N, ...
                                0, ...
                                @InventoryStateIterator, ...
                                @InventoryDecisionIterator, ...
                                @InventoryTerminalValueIterator, ...
                                @InventoryCostToGoFn, ...
                                @InventoryStateTransistionProbability, ...
                                parameters);

toc

%% ========== A small test case ==========

fprintf('\n\n---- Small ADP test case ----\n')

% Test Parameters
max_inv = 5;
n_periods = 4;
lambda = 1;
disc_rate = 0.05;

adp_samples = 200;

% First run exact dp for comparison
disp('Run DP basecase for comparision')
tic
[small_DP_orders, small_DP_values] = inventory3(max_inv, n_periods, lambda, disc_rate);
toc

%Now run adp version
fprintf('\nRun ADP version with %d samples\n', adp_samples)
tic
[small_ADP_orders, small_ADP_results, small_ADP_values] = inventory_adp(n_periods, adp_samples, ...
    struct('lambda', lambda, 'max_inv', max_inv, 'dr', disc_rate), ... %inv_problem structure
    struct('plot', 0, 'bootstrap', 3, 'smooth', [2 2], 'trace', true));  %adp_setup structure
toc
fprintf('\nDifference of optimal inventory state:')
small_ADP_orders - small_DP_orders %#ok<NOPTS,MNEFF>


%% ========== Medium problem ==========
%Note: similar to the simple call above, but using fancier options

%params
N = 5;  %Number of periods
iter = 1000;

fprintf('\n\n---- Medium ADP test case ----\n')

% First run exact dp for comparison
disp('DP basecase for comparision')
fprintf('  n_periods = %d\n', N)
tic
[med_DP_orders, med_DP_values] = inventory3(50, N, 20, 0.05);
toc


fprintf('\nRun ADP version with %d samples\n', iter)
fprintf('  n_periods = %d, iterations = %d\n', N, iter)

tic
[med_ADP_orders, med_ADP_results, med_ADP_values] = ...
    inventory_adp(N, iter, [], ...
        struct('bootstrap', N/10, 'trace',true, 'plot', false, 'fix_rand', true ...
               ,'stepsize', '1overN' ...
               ...,'step_opt', 10 ...
               ...,'stepsize', @ssSTC,'step_opt', struct('step0', 1, 'a', 1, 'b', 1, 'beta', 1)...
        ),...
        [],med_DP_orders);
toc

%Compute and display errors
figure
subplot(2,1,1)
plot(med_ADP_results.trace.rel_err)
title('Relative Error')
subplot(2,1,2)
plot(med_ADP_results.trace.abs_err)
title('Absolute Error')

med_ADP_order_error = med_DP_orders - med_ADP_orders;
fprintf('   %g%% Error in decision space (ratio Frobenius norm (RMS) errors:optimal)\n', ...
             norm(med_ADP_order_error, 'fro')/norm(med_DP_orders, 'fro') * 100)
fprintf('   RMS Error by period:')
for period = 1:size(med_ADP_order_error, 2);
    fprintf(' %g', norm(med_ADP_order_error(:,period)))
end

fprintf('\n')
fprintf('   Absolute decision error max: %d (%%%g)\n', ...
    max(max(abs(med_ADP_order_error))), max(max(abs(med_ADP_order_error)))/max(max(med_DP_orders))*100)
fprintf('       by period (max opt order = %d):',max(max(med_DP_orders)))
fprintf(' %d', max(abs(med_ADP_order_error)))
fprintf('\n')     
%% ========== A larger comparision ==========
% Note: this one needs some tuning still
% This also demonstrates the use of Inv_Prob struct compatability mode for
% inventory_dp

% Build struct of Test Parameters
big_inv_prob.max_inv = 15;
big_inv_prob.n_periods = 10;
big_inv_prob.lambda = ceil(big_inv_prob.max_inv /2);
big_inv_prob.dr = 0.05;
big_inv_prob.order_cost_fix = 4;

adp_samples = 2000;      %Note: this will be max samples if tol specified
num_adp_runs = 5;

fprintf('\n\n---- Big ADP test case ----\n')
fprintf('  max_inv = %d, lambda = %g, n_periods = %d\n', ...
    big_inv_prob.max_inv, big_inv_prob.lambda, big_inv_prob.n_periods)

% First run exact dp for comparison but only if 0 or 1 ADP run is planned
if not(exist('big_DP_orders', 'var'))  ...
        || not(isequal(size(big_DP_orders), [big_inv_prob.max_inv + 1, big_inv_prob.n_periods])) ...
        || num_adp_runs <= 1
    disp('Run DP basecase for comparision')
    tic
    [big_DP_orders, big_DP_values] = inventory_dp(big_inv_prob);
    toc
end

%adp_setup structure:
% Note: inv4err set to compute error only for the 95% likely pre-decision
% states based on the DP optimal post decsion state - 2.5 percental of the
% poisson cdf. (poissinv(.025, lambda))
big_inv_adp = struct(   'plot',         0, ...
                    	'bootstrap',    max(min(2,adp_samples/10),50)  ...
                        ,'smooth',      [5 5] ...
                        ...,'trace',       true...  %enabled by default when tol specified
                        ,'inv4err',     0:ceil(big_DP_orders(1,1)-poissinv(.025, big_inv_prob.lambda)) ...
                        ...,'inv4err',     0:32 ...    % 95% likely pre decision states for max_inv = 250, lambda = 125
                        ...,'inv4err',     0:46 ...    % 95% likely pre decision states for max_inv = 500, lambda = 250
                        ...,'inv4err',     0:65 ...    % 95% likely pre decision states for max_inv = 1000, lambda = 500
                        ...,'inv4err',     0:93 ...    % 95% likely pre decision states for max_inv = 2000, lambda = 1000
                        ...,'fix_rand',    true ...
                        ...,'stepsize',    '1overN' ...   %Default=harmonic (can be much faster)
                        ...,'stepsize',    'STC2C' ...   %Default=harmonic (can be much faster)
                        ...,'step_opt',    struct('a', 2, 'b', 1000, 'target', 0.0005, 'step2', 1, 'beta', 1)...          %step size algorithm parameter(s)
                      );
big_tol = 0.02;     %Policy convergence toleranace (Frobenius norm sense)
    
%Initialize ADP sample results storage
big_ADP_time = zeros(1,num_adp_runs);
big_ADP_iter = zeros(1,num_adp_runs);
big_ADP_err = zeros(1,num_adp_runs);

%Now run adp version
fprintf('\nRun ADP version with max %d passes\n', adp_samples)

for r = 1:num_adp_runs
    fprintf('\nADP trial #%d...\n', r)
    tic
    [big_ADP_orders, big_ADP_results, big_ADP_values] ...
        = inventory_adp(big_inv_prob.n_periods, adp_samples, big_inv_prob,... 
                        big_inv_adp, [], big_DP_orders, big_tol);
    toc

    %Store results
    big_ADP_time(r) = big_ADP_results.trace.time_wo_trace(big_ADP_results.n);
    big_ADP_iter(r) = big_ADP_results.n;
    big_ADP_err(r) = big_ADP_results.trace.abs_err(big_ADP_results.n);
    
    
    %Compute and display errors
    if isfield(big_ADP_results, 'trace')
        fprintf('ADP run time, not counting trace is %g seconds\n', ...
            big_ADP_time(r))

%         figure(r)
%         subplot(2,1,1)
%         plot(big_ADP_results.trace.rel_err)
%         title(sprintf('Relative Error (inv=%d)',big_inv_prob.max_inv))
%         subplot(2,1,2)
%         plot(big_ADP_results.trace.abs_err)
%         title('Absolute Error')
    end

    big_ADP_order_error = big_DP_orders(big_inv_adp.inv4err+1,:) - big_ADP_orders(big_inv_adp.inv4err+1,:);
    fprintf('   %g%% Error in likely decision space (ratio Frobenius norm (RMS) errors:optimal)\n', ...
                 norm(big_ADP_order_error, 'fro')/norm(big_DP_orders, 'fro') * 100)
    fprintf('   RMS Error by period:')
    for period = 1:size(big_ADP_order_error, 2);
        fprintf(' %g', norm(big_ADP_order_error(:,period)))
    end

    fprintf('\n')
    fprintf('   Absolute decision error max: %d (%g%%)\n', ...
        max(max(abs(big_ADP_order_error))), max(max(abs(big_ADP_order_error)))/max(max(big_DP_orders))*100)
    fprintf('       by period (max opt order = %d):',max(max(big_DP_orders)))
    fprintf(' %d', max(abs(big_ADP_order_error)))
    fprintf('\n')
end

fprintf('\nAverage of %d ADP runs: \n', num_adp_runs)
fprintf(  '----------------------- \n')
fprintf(  '   Time (wo trace): %g sec\n', mean(big_ADP_time))
fprintf(  '   Num iterations: %g\n', mean(big_ADP_iter))
fprintf(  '   Abs (policy) error: %g%%\n', mean(big_ADP_err)*100)
