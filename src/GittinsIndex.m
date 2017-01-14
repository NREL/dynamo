function [ g_idx ] = GittinsIndex( num_obs, disc_factor, var_is_known, means, stdevs )
%GITTINSINDEX Gittins Index for normally distributed rewards in multi-armed bandit problems
%
%   g_idx = GittinsIndex(num_obs, disc_rate, var_is_known, means, stdevs)
%
% INPUTS
%  num_obs       vector of number of observations of a bandit process
%  disc_factor   scalar discount factor (default = 0.95) IMPORTANT: only
%                0.95 and 0.99 implemented
%  var_is_known  boolean scalar. true if variance is known (default=true)
%  means         vector of means corresponding to num_obs (default = 0)
%  stdevs        vector of standard deviations corresponding to num_obs
%                (default = 1)
%
% Notes:
%  - Use [] for early parameters to keep defaults when specifying later
%    options
%
% IMPORTANT:
%  - This version currently relies on the look-up table from Powell's book
%    (2007, p338, table 10.2) As such it is limited to the values in that
%    table, most notably only 0.95 and 0.99 discount factors are supported
%  - indicies for num_obs > 100 values decrease from the table value of 100
%    by a factor of 1/sqrt(n)
%  - Obs = 0 (and 1 for unknown variance) return Inf gittins indicies
%
% Reference:
%  Powell, W. B. (2007). Exploration vs Exploitation. In Approximate
%    Dynamic Programming: Solving the Curses of Dimensionality
%    (1st ed., pp. 323-350). Hoboken, NJ (US): Wiley-Interscience.
%
% MATLAB implementation originally by Bryan Palmintier, 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-10-18 15:14  BryanP      Initial version
%   2  2010-10-18 20:30  BryanP      Added support for empty early parameters

%-- Handle options
if nargin < 2 || isempty(disc_factor)
    disc_factor = 0.95;
end

if (disc_factor ~= 0.99) && (disc_factor ~= 0.95)
    error('Only discount factors of 0.99 and 0.95 are supported')
end

if nargin < 3 || isempty(var_is_known)
    var_is_known = true;
end

if nargin < 4 || isempty(means)
    means = zeros(size(num_obs));
end

if nargin < 5 || isempty(stdevs)
    stdevs = ones(size(num_obs));
end

%-- Setup Lookup table values
% Important: obs must be in order and start with 1
% Powell (2007) Table 10.2, p 338
%               -- known var -  - unknown var -
%         obs    0.95    0.99    0.95     0.99
table = [  1    0.9956	1.5758     Inf     Inf
           2    0.6343  1.0415 10.1410 39.3343
           3    0.4781  0.8061  1.1656  3.1020
           4    0.3878  0.6677  0.6193  1.3428
           5    0.3281  0.5747  0.4478  0.9052
           6    0.2853  0.5072  0.3590  0.7054
           7    0.2528  0.4554  0.3035  0.5901
           8    0.2274  0.4144  0.2645  0.5123
           9    0.2069  0.3808  0.2353  0.4556
          10    0.1889  0.3528  0.2123  0.4119
          20    0.1058  0.2094  0.1109  0.2230
          30    0.0739  0.1520  0.0761  0.1579
          40    0.0570  0.1202  0.0582  0.1235
          50    0.0464  0.0998  0.0472  0.1019
          60    0.0392  0.0855  0.0397  0.0870
          70    0.0339  0.0749  0.0343  0.0760
          80    0.0299  0.0667  0.0302  0.0675
          90    0.0267  0.0602  0.0269  0.0608
         100    0.0242  0.0549  0.0244  0.0554
        ];
obs = table(:,1);
if var_is_known
    if disc_factor == 0.95
        g_idx_lookup = table(:,2);
    elseif disc_factor == 0.99
        g_idx_lookup = table(:,3);
    end
else
    if disc_factor == 0.95
        g_idx_lookup = table(:,4);
    elseif disc_factor == 0.99
        g_idx_lookup = table(:,5);
    end
end

%-- Compute the indicies
% First find which values are outside the range of the table (index vector)
extrap = num_obs > obs(end);

% lookup the standard normal indicies that are within the table
g_idx(not(extrap)) = interp1(obs, g_idx_lookup, num_obs(not(extrap)));

% fix any NaNs introduced by interp1 not liking infinity for the obs=1
g_idx(num_obs == 1) = g_idx_lookup(1);

% fix any NaNs introduced by interp1 for zero obs by setting to Inf
g_idx(num_obs == 0) = Inf;

% extrapolate higher values of n using 1/sqrt(n)
g_idx(extrap) = g_idx_lookup(end) .* sqrt(obs(end)) ./ sqrt(num_obs(extrap));

% then scale to our mean and std deviation
g_idx = means + stdevs .* g_idx;

end
