function compareStepSize(n_max, a, target, b, step2, half_jump_step)
%COMPARESTEPSIZE Test plots to compare stepsize formulas
%
% Usage:
%   results = compareStepSize(n_max, a, target, b, step2, half_jump_step)

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-03-15 13:40  BryanP      Initial Code
%   2  2012-01-03 15:25  BryanP      Added step2 parameter 

if nargin < 5 || isempty(step2)
    step2 = 1;
end

if nargin <6
    half_jump_step = [];
end

%Setup
names = {
         'McClain'
         'Harmonic2C'
         '1overN'
         'Harmonic'
         'STC'
         'STC2C'
         '1overNrecur'
         'HarmonicRecur'
         'STCrecur'
         'STC2Crecur'};
     
n_names = length(names);
recur_versions = 4;  %number of final formulas that are recursive duplicates and should be plotted using contrasting dotted lines.

%Initialize
results = ones(length(names), n_max);

func = cell(size(names));
for f = 1:n_names
    func{f} = str2func(['ss' names{f}]);
end

%-- Setup the remaining functions that all take a param structure
param.a = a;
param.target = target;
param.b = b;
param.step2 = step2;
param.beta = 1;

for f = 1:n_names
    results(f, 1) = func{f}(1, [], param);
    param.step2 = step2;
    for n = 2:n_max;
        if isempty(half_jump_step) || n ~= floor(n_max/2)
            old_step = results(f, n-1);
        else
            old_step = half_jump_step;
            param.step2 = half_jump_step;
        end
        results(f, n) = func{f}(n, old_step, param);
    end
end
 
%-- Now make the plot
%clear existing figure
clf
%plot main formulas
plot(1:n_max, results(1:(end-recur_versions),:), 'LineWidth',1.1)
%overplot recursive versions as dashed lines. Initial formulas that have no
%recursive equivalent ensure a collor shift
hold on
plot(1:n_max, results((end-recur_versions+1):end,:),':', 'LineWidth', 1.1)
legend(names)
%Build up title string
title_str = sprintf('Stepsize comparison for a=%d, b=%d, target=%g',a,b,target);
title_str = [title_str sprintf(', 2nd step=%g', step2)];
if not(isempty(half_jump_step))
    title_str = [title_str sprintf(', step jump=%g', half_jump_step)];
end

title(title_str)
xlabel('n')
ylabel('stepsize')