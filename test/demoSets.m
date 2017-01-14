%Simple test/demo for Sets

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   2  2016-09-29 14:30  BryanP      Updated to match setBasic v3 (and v4) 
%   1  2016-07-07        BryanP      Initial version 



%% Basic sets
test = setBasic([0 1.5 3 3.2; 2 3 4.2 4], 'sobol', [.5 .3 .2 0])
s_list = test.as_array();  % Get full list of states (all combinations): Big!
size(s_list)            % Should be 2310 x 4
s_list(1:10,:)          % First 10 states

test.SampleType         % Should be 'sobol'
test.sample()           % Single sample
test.sample(10)         % Ten samples
disp('Next set of samples will be psuedorandom, not sobol')
test.sample(10,'rand')  % Should switch to psuedo-random numbers
test.SampleType         % Should be 'rand'
test.sample_init('sobol')
test.SampleType         % Should be 'sobol'
test.sample()           % Single sample (Now using sobol... Should be the same as earlier samples)
test.sample(10,'sobol') % Should throw a warning and still return the samples (same as earlier sobol)
test.SampleType         % Should be 'sobol'
test.sample(10)         % Ten samples

range = test.range()            % Range of values for all dimensions (matches input)
