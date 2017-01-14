% Scratchpad for:
% modularized Dynamic Programming implementation of MULTI PRODUCT Inventory problem
% Developed as part of ESD.862 (MIT) Final Project
% Bryan Palmintier
% Spring 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
% many 2010              BryanP      original code, additions not tracked
%   2  2010-05-26 23:50  BryanP      MultiInvCombin tests
%   3  2010-08-07 11:50  BryanP      split from InvDP_scratchpad
%   4  2011-03-29 10:10  BryanP      Renamed MultiInvValidCombin to CombinWithLimits
%   5  2011-05-06 15:40  BryanP      Added Parallel DP comparison

% Adapted from hw4_scratchpad.m, which can also be used to check results.
% Can also compare results with InvDP_scratchpad

%% ------ Multi-Product DP Testing ------

%% --- Test constrained state space generation
%Illustrate some of the features of the constrained combination tool:
                                                            % size min max
%A large 3 item space only constrained by max per item (10)
% Note: zero inventory state OK, so have 11x11x11
a=CombinWithLimits([1 1 1], Inf, 10);                    % 1331   0  10
%Per item max inventory, and default for space constraint
b=CombinWithLimits([1 1 1], [], [8 9 10]);               %  990   0  10
%Add total space (warehouse size) constraint.
c=CombinWithLimits([1 1 1], 15, [8 9 10]);               %  641   0  10
%Consider products with non-uniform sizes
d=CombinWithLimits([1 2 3], 15, [8 9 10]);               %  151   0   8
%Add a minimum total inventory constraint
e=CombinWithLimits([1 2 3], 15, [8 9 10], 14);           %   39   0   8
%Include minimum inventory per item
f=CombinWithLimits([1 2 3], 15, [8 9 10], 14, [3,2,1]);  %    9   1   8

%--A few examples with backorders...
%Setting a neg doesn't work b/c default to no backorders (same list as d)
g=CombinWithLimits([1 2 3], 15, [8 9 10], -5);           %  151   0   8
%Have to explicitly set negative total and item quantities
h=CombinWithLimits([1 2 3], 15, [8 9 10], -5,-5);        %  987  -5   8
%and can have max per product backorders
i2=CombinWithLimits([1 2 3], 15, [8 9 10], -5,-2);       %  479  -2   8
%which saves more space than having a smaller max net backorder list 
j2=CombinWithLimits([1 2 3], 15, [8 9 10], -2,-5);       %  835  -5   8
%We can also create entirely negative lists
k=CombinWithLimits([1 2 3], -5, [], -20,-5);             %    3  -5  -2

%% --- Test reverse lookup performance
% Note: change commenting in MultiInvLookup() to compare algorithms

params = MultiInvInit(50,[7 3 1 4],[1 2 3 4]);

%test = randi(1e5,1,[1 params.n_states]);test_s=params.states(test,:);
test = randi(params.n_states,[params.n_states,1]);test_s=params.states(test,:);
tic,idx=MultiInvLookup(test_s,params);toc,nnz(idx ~= test)

test = (1:2:358)';test_s=params.states(test,:);
tic,idx=MultiInvLookup(test_s,params);toc,nnz(idx ~= test)

%% --- Treat the super simple single product case as multi-product as a check
% Note: the results for the first 5 states should match DiscDP_Orders{2}
% and DiscDP_Values{2}... and they do!
p_demand = [0.25, 0.5, 0.25];
max_inv = 3;
time_steps = 3;
disc_rate = 0;

simplemulti_params=MultiInvInit(max_inv,[1 1],{1;p_demand});

fprintf('\nStarting MultiInv framework for a single product...')
tic
[SimpleMulti_Orders, SimpleMulti_Values] = DP(time_steps, disc_rate, ...
    @MultiInvState, @MultiInvDecision, @MultiInvTermValue, @MultiInvCost, ...
    @MultiInvTransProb, simplemulti_params);
toc

%% --- Now check a bit larger single product as a fake multi product...
% Note: the results for the first max_inv states should match DiscDP_Orders{2}
% and DiscDP_Values{2}... and they do!
max_inv = 10;
lambda = 5;
time_steps = 15;
disc_rate = 0.1;

testmulti_params=MultiInvInit(max_inv,[1 1],[0 lambda]);

fprintf('\nStarting dummy multi-product inventory with zero demand for one product...')
tic
[TestMulti_Orders, TestMulti_Values] = DP(time_steps, disc_rate, ...
    @MultiInvState, @MultiInvDecision, @MultiInvTermValue, @MultiInvCost, ...
    @MultiInvTransProb, testmulti_params);
toc

%% --- Run a medium-sized 3 product case
% Note this is a multi product inventory model with storage space
% constraints and independent poisson demands.
multi3_params=MultiInvInit(20,[2 3 1],[1 2 3]);
time_steps = 5;
disc_rate = 0.1;

fprintf('\nStarting medium-sized 3 product case (may take over a minute)...\n')
tic
[Multi3_Orders, Multi3_Values] = DP(time_steps, disc_rate, ...
    @MultiInvState, @MultiInvDecision, @MultiInvTermValue, @MultiInvCost, ...
    @MultiInvTransProb, multi3_params,50);
toc

%% Compare to parallel DP version
% fprintf('\nParallel Programming for same medium-sized 3 product case...\n')
% 
% if matlabpool('size') == 0
%     %Setup to use the default (maximum) number of cores
%     fprintf('  Enabling parallel workers...\n')
%     tic
%     matlabpool open
%     toc
% else
%     n_cores = matlabpool('size');
%     fprintf('  Using %d previously initialized parallel workers\n', n_cores)
% end

tic
[ParMulti3_Orders, ParMulti3_Values] = parDP(time_steps, disc_rate, ...
    @MultiInvState, @MultiInvDecision, @MultiInvTermValue, @MultiInvCost, ...
    @MultiInvTransProb, multi3_params,50);
toc
