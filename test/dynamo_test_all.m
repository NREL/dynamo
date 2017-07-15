function dynamo_test_all
% Run all tests for dynamo
%  Add additional tests, doctest or otherwise, as needed
%
%  IMPORTANT: MATLAB changed their output formatting in R2016b for all
%  non-double array displays. as a result many of these doctests will not
%  work. 
% 
%  Note: `example/multi_inventory` needs to be in the path

%% ====== DOCTESTS =======
% Warn if using an older MATLAB than R2016b. Note version 9.1 = R2016b
if verLessThan('matlab','9.1')
    warning('dynamo_run_all_tests:matlab_version_too_old', ...
        'MATLAB changed its output formatting in R2016b. Many of these doctests will break on older versions')
end

test_time = tic;

docs_to_test = {    
                    %Sets
                    'setSingleItem'
                    'setCombinWithLimits'
                    %Random Processes
                    'rpDiscreteSample'
                    'rpLattice'
                    %Utilities
                    'utilRandSetCurState'
                    'utilRandStatefromState'
                    %Examples
                    'doctest_MultiInv'
                    'MultiInv_demo'
                    'MultiInvSetupProblem'
                };
%initialize pass & test counts
total_tests = 0;
total_pass = 0;

for d_idx = 1:length(docs_to_test)
    d_name = docs_to_test{d_idx};
    
    fprintf('\nTesting %s\n', d_name)
    [~, n_pass, n_tests] = doctest(d_name);
    total_tests = total_tests + n_tests;
    total_pass = total_pass + n_pass;
end

%Display summary
fprintf('\n  DOCTEST Results (%d files): ', length(docs_to_test))
if total_pass == total_tests
    fprintf('PASS')
else
    fprintf('FAIL')
end
fprintf(', %d/%d tests pass (%d%%)\n  ', total_pass, total_tests, total_pass/total_tests*100)
toc(test_time)

%% === MATLAB Unit Tests ===
test_results = runtests('testrpLattice.m');
