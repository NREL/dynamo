% Cavetest: Working development script for the CAVE algorithm (ESD.862)
%
% Reference: Godfrey & Powell, 2001
%

% HISTORY
% ver     date    time        who      changes made
% ---  ---------- -----  ------------- ---------------------------------------
%   1  2010-04-18 01:14  BryanP        Split from CaveUpdate ver 6
%   2  2010-04-18 03:00  BryanP        Added parabola test case
%   3  2010-04-20 19:00  BryanP        Streamlined using CaveTest
%   4  2010-04-21 15:30  BryanP        Added quick non-convexity checks

%% --- Example #1 Natural Log with flat start & 20 samples
figure(1)
clf

disp('Natural Log with Flat Start')
[u,v,s]=CaveTest(@log,[1 40],20);

% display some quick results & quick check
samples=s'
nonconvexities=nnz(find(diff(v)>0))

figure(1)

%% --- Example 2: Natural Log with Mort's simple estimate & n samples

disp('Natural Log with Mort''Start')
%parameters
x_max = 40;
x_min = 1;
ds = 0.1;
my_slopes = @(s) [(log(s)-log(s-ds))/ds, (log(s+ds)-log(s))/ds];

n_samples = 20;

%plot actual values
x=x_min:x_max;
y=log(x);

figure(2)
clf
plot(x,y,'g','LineWidth',2)
title(sprintf('natural log with Mort''s start & %d random samples',n_samples))

% Initial Approximation
u=[   1   10   20   30 40];
v=[0.35 0.25 0.15 0.01 0];

% Samples
obs = rand(n_samples,1)*(x_max-x_min) + x_min

%Actually run the algorithm & make pretty pictures
[u,v] = CaveUpdate(obs, my_slopes, u, v, [x_min x_max], 1, [2 2], [], 1);

% display some quick results & quick check
nonconvexities=nnz(find(diff(v)>0))
figure(2)

%% --- Ex3: Our sample case from the old version of cavetest 
disp('Natural Log with Mort''s Start and Sarvee''s samples (scrambled)')
%parameters
x_max = 40;
x_min = 1;
ds = 0.1;
my_slopes = @(s) [(log(s)-log(s-ds))/ds, (log(s+ds)-log(s))/ds];

n_samples = 20;

%plot actual values
x=x_min:x_max;
y=log(x);

figure(3)
clf
plot(x,y,'g','LineWidth',2)
title(sprintf('natural log with Mort''s start & Sarvee''s samples',n_samples))
% Initial Approximation
u=[   1   10   20   30 40];
v=[0.35 0.25 0.15 0.01 0];

obs = [35     6    25     3    16];
%scramble observations
obs=obs(randperm(length(obs)))

%Actually run the algorithm & make pretty pictures
[u,v] = CaveUpdate(obs, my_slopes, u, v, [x_min x_max], 1, [2 2], [], 1);

% display some quick results & quick check
nonconvexities=nnz(find(diff(v)>0))
figure(3)

%% --- Ex4: An inverse parabola y = -(x/10-2)^2 + 4 with flat start & n samples
disp('Inverted Parabola with Flat Start')
%parameters
n_samples = 100;
my_fun = @(x) -(x./10-2).^2 + 4;

figure(4)
clf
[u,v,s]=CaveTest(my_fun,[0 40],n_samples,'inverted parabola');

% display some quick results & quick check
samples=s'
nonconvexities=nnz(find(diff(v)>0))
figure(4)