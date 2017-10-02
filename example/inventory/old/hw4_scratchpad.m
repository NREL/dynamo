% HW 4 Scratchpad for ESD.862
% Theme: Dynamic Programming
% Bryan Palmintier
% Spring 2010

%% ========== Question 1 - Inventory Problem (Finite Horizon) ==========

%% --- Part 1: Repeat Putterman 3.2.2
[Simple_Orders, Simple_Values] = inventory_dp(3, 3, [0.25, 0.5, 0.25]);

%% --- Part 2: Larger system with various Poisson demands
time_steps = 15;
max_inv = 10;
run_plots = false;
max_lambda = 20;

fprintf('\nRunning custom inventory DP for %d lambda values\n',max_lambda)
tic
for lambda = 1:max_lambda;
	[Finite_Orders{lambda}, Finite_Values{lambda}] = inventory_dp(max_inv, time_steps, lambda);

    if run_plots
        %Plot Values
        figure(lambda*2-1)
        surf(Finite_Values{lambda})
        xlim([1,time_steps+1])
        xlabel('Time')
        ylim([1, max_inv+1])
        ylabel('Inventory State')
        title(sprintf('Values of State over time for lambda = %d',  lambda))

        %Plot Orders
        figure(lambda*2)
        surf(Finite_Orders{lambda})
        xlim([1,time_steps+1])
        xlabel('Time')
        ylim([1, max_inv])
        ylabel('Inventory State')
        title(sprintf('Optimal Order by state and time for lambda = %d',  lambda))
    end
end
toc

%% --- Part 3: Add Discounting
time_steps = 15;
max_inv = 10;
lambda = 5;

disc_rates = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5];

clear Disc_Orders Disc_Values Disc_p_demand
for dr = 1:length(disc_rates);
	[Disc_Orders{dr}, Disc_Values{dr}, Disc_p_demand{dr}] ...
        = inventory_dp(max_inv, time_steps, lambda, disc_rates(dr));
end

%% ========== Question 2 - Infinite Horizon Inventory Problem ==========

%% --- Part 1: Value Iteration AND Part 2: Policy Iteration
lambdas = [   1    3    6    5    5    5    5    5    5];
drs     = [0.05 0.05 0.05 0.05  0.1  0.2 0.05 0.05 0.05];
tols    = [  -4   -4   -4   -4   -4   -4    0    -2  -6];
tols = 10.^tols;

clear Val_Iter_Orders Val_Iter_Values Val_Iter_n_iter
for iter = 1:length(lambdas)
    fprintf(1,'\rRunning iteration #%d:',iter)
    fprintf(1,'Value iteration...')
    [Val_Iter_Orders(:,iter), Val_Iter_Values(:,iter), Val_Iter_n_iter(:,iter)] = ...
        inventory_inf_horiz(10, lambdas(iter), drs(iter), 'value', tols(iter));
    fprintf(1,'Policy iteration...')
    [Pol_Iter_Orders(:,iter), Pol_Iter_Values(:,iter), Pol_Iter_n_iter(:,iter)] = ...
        inventory_inf_horiz(10, lambdas(iter), drs(iter), 'policy');
end
disp('Done')

%% --- Part 2b: Trace Results
[Orders, Values, n_iter, Val_Iter_Trace] = inventory_inf_horiz(10, 5, 0.05, 'value',10);
[Orders, Values, n_iter, Pol_Iter_Trace] = inventory_inf_horiz(10, 5, 0.05, 'policy');

%pickout selected value function (recall that the state is one more than
%the inventory level
trace_state = 3;

figure(13)
clf
hold on
plot(1:size(Val_Iter_Trace,2),Val_Iter_Trace(trace_state,:),...
    1:size(Pol_Iter_Trace,2),Pol_Iter_Trace(trace_state,:), ...
    'LineWidth',2)
plot([0 size(Val_Iter_Trace,2)],Pol_Iter_Trace(trace_state,end) * [1 1],'k--')
title(sprintf('Convergence results for inventory = %d value function', trace_state-1))
xlabel('iteration number')
legend('Value Iteration', 'Policy Iteration', 'Location','SouthEast')