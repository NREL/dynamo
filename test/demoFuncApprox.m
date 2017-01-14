%Quick & dirty demonstration of different FuncApprox
%
% TODO: 
%   -- add faLocalAvg
%   -- use different function: peaks? 

clc
clear

%% Create some fake data
%initial values
pts = bsxfun(@times, randi(10,10,2), [1 1]);
vals = rand(10,1);

% values to mix-in, with some duplicate points
% new_pts = vertcat(pts(1:3,:),rand(7,2));
% new_vals = rand(10,1);

c = faDiscrete('pts', pts, 'vals', vals, 'max_size_per_dim', [10, 10]);

% c.plot

%b = faInterp(pts, vals), figure(2), b.plot, title('(natural) Interpolation (n=10)');

% %% Create and plot function approximations based on initial data points
% %a = faThinPlate(pts, vals), figure(1), a.plot, title('Thin Plate (n=10)')
% b = faInterp(pts, vals), figure(2), b.plot, title('(natural) Interpolation (n=10)')
% c = faDiscrete(pts, vals);
% 
% %Use 3 nearest, b/c default ball of 0.25 usually doesn't have enough points
% k = faLocalRegr(pts, vals), figure(3), k.plot(3), title('Local Regression -- 3 nearest (n=10)')
% 
% %% Add in new values & plot again
% a.update(new_pts, new_vals), figure(4), a.plot, title('Thin Plate (n=20)')
% b.update(new_pts, new_vals), figure(5), b.plot, title('(natural) Interpolation (n=20)')
% k.update(new_pts, new_vals), figure(6), k.plot(3), title('Local Regression -- 3 nearest  (n=20)')
% 
% %% Alternative approximations
% figure, b.plot('nearest'), title('(nearest) Interpolation (n=20)')
% figure, b.plot('linear'), title('(linear) Interpolation (n=20)')
% 
% figure, k.plot(1-eps, true), title('Local Regression -- all, dist_weight  (n=20)')
% figure, k.plot(1), title('Local Regression -- only 1 nearest  (n=20)')
% figure, k.plot(2), title('Local Regression -- 2 nearest  (n=20)')
% figure, k.plot(5), title('Local Regression -- 5 nearest  (n=20)')
% figure, k.plot(5, true), title('Local Regression -- 5 nearest, dist_weight  (n=20)')
% 
new_pts = vertcat(pts(1:3,:),rand(7,2));
new_vals = rand(10,1);

%% Create and plot function approximations based on initial data points
a = faThinPlate(pts, vals), figure(1), a.plot, title('Thin Plate (n=10)')
b = faInterp(pts, vals), figure(2), b.plot, title('(natural) Interpolation (n=10)')

%Use 3 nearest, b/c default ball of 0.25 usually doesn't have enough points
k = faLocalRegr(pts, vals), figure(3), k.plot(3), title('Local Regression -- 3 nearest (n=10)')

%% Add in new values & plot again
a.update(new_pts, new_vals), figure(4), a.plot, title('Thin Plate (n=20)')
b.update(new_pts, new_vals), figure(5), b.plot, title('(natural) Interpolation (n=20)')
k.update(new_pts, new_vals), figure(6), k.plot(3), title('Local Regression -- 3 nearest  (n=20)')

%% Alternative approximations
figure, b.plot('nearest'), title('(nearest) Interpolation (n=20)')
figure, b.plot('linear'), title('(linear) Interpolation (n=20)')

figure, k.plot(1-eps, true), title('Local Regression -- all, dist_weight  (n=20)')
figure, k.plot(1), title('Local Regression -- only 1 nearest  (n=20)')
figure, k.plot(2), title('Local Regression -- 2 nearest  (n=20)')
figure, k.plot(5), title('Local Regression -- 5 nearest  (n=20)')
figure, k.plot(5, true), title('Local Regression -- 5 nearest, dist_weight  (n=20)')

