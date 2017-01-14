% Scratchpad for:
% modularized Dynamic Programming implementation of SINGLE PRODUCT Inventory
% Developed as part of ESD.862 Final Project
% Bryan Palmintier
% Spring 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
% many 2010              BryanP      original code, additions not tracked
%   2  2010-05-26 23:50  BryanP      MultiInvCombin tests
%   3  2010-08-07 11:50  BryanP      Removed MultiInv to focus only on Inv

% Requires these functions:
% DP, InvState, InvDecision, InvTermValue, InvCost, InvTransProb

% Adapted from hw4_scratchpad.m, which can also be used to check results

%% ========== Question 1 - Inventory Problem (Finite Horizon) ==========

%% --- Part 1: Repeat Putterman 3.2.2
p_demand = [0.25, 0.5, 0.25]';
max_inv = 3;
%get problem defaults & compute poisson demand distribution
simple_params = InvInit(max_inv,p_demand);

[SimpleDP_Orders, SimpleDP_Values] = DP(3, 0, @InvState, @InvDecision, ...
    @InvTermValue, @InvCost, @InvTransProb, simple_params);

%% --- Part 2: Larger system with various Poisson demands
time_steps = 15;
max_lambda = 20;
max_inv = 10;

fprintf('\nRunning modularized DP for inventory problem with %d lambda values\n', max_lambda)
clear FiniteDP_Orders FiniteDP_Values
tic
for lambda = 1:max_lambda;
    %get problem defaults & compute poisson demand distribution
    finite_params = InvInit(max_inv, lambda);

    [FiniteDP_Orders{lambda}, FiniteDP_Values{lambda}] = DP(time_steps, 0, ...
        @InvState, @InvDecision, @InvTermValue, @InvCost, @InvTransProb, ...
        finite_params); %#ok<SAGROW>
end
toc

%% --- Part 3: Add Discounting
time_steps = 15;
max_inv = 10;
lambda = 5;
%get problem defaults & compute poisson demand distribution
disc_params = InvInit(max_inv, lambda);

disc_rates = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5];

clear Disc_Orders Disc_Values
for dr = 1:length(disc_rates);
	[DiscDP_Orders{dr}, DiscDP_Values{dr}] = DP(time_steps, disc_rates(dr), ...
        @InvState, @InvDecision, @InvTermValue, @InvCost, @InvTransProb, ...
        disc_params); %#ok<SAGROW>
end 

