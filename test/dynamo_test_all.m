function dynamo_test_all
% Run all tests for dynamo
%  Add additional tests, doctest or otherwise, as needed
%
%  IMPORTANT: MATLAB changed their output formatting in R2016b for all
%  non-double array displays. as a result many of these doctests will not
%  work. 
% 
%  Note: `example/multi_inventory` needs to be in the path

% Warn if using a newer MATLAB than R2016a. Note version 9.1 = R2016b
if not(verLessThan('matlab','9.1'))
    warning('dynamo_run_all_tests:matlab_version_too_new', ...
        'MATLAB changed their output formatting in R2016b. Many of these doctests will break')
end

docs_to_test = {    
                    %Sets
                    'setSingleItem'
                    'setCombinWithLimits'
                    %Random Processes
                    'rpDiscreteSample'
                    %Utilities
                    'utilRandSetCurState'
                    'utilRandStatefromState'
                    %Examples
                    'doctest_MultiInv'
                    'MultiInv_demo'
                    'MultiInvSetupProblem'
                };

for d_idx = 1:length(docs_to_test)
    d_name = docs_to_test{d_idx};
    
    fprintf('Testing %s\n', d_name)
    doctest(d_name)
end

