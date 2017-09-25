function dynamo_test_all
% Run all tests for dynamo
%  Add additional tests, doctest or otherwise, as needed
%
%  IMPORTANT: MATLAB changed their output formatting in R2016b for all
%  non-double array displays. as a result many of these doctests will not
%  work. 
% 
%  Note: `example/multi_inventory` needs to be in the path

%% === MATLAB Unit Tests ===
fprintf('\n------ UNIT TESTS -------\n')
unittest_list = {    
                    %Sets

                    %Random Processes
                    'testrpLattice.m'
                    %Utilities

                    %Examples

                    };
%initialize pass & test counts
total_unittests = 0;
total_unitpass = 0;

unittest_time = tic;
for u_idx = 1:length(unittest_list)
    u_name = unittest_list{u_idx};
    
    test_results = table(runtests('testrpLattice.m'));
    n_pass = sum(test_results.Passed);
    n_tests = n_pass + sum(test_results.Incomplete) + sum(test_results.Failed);
    
    fprintf('%s: ', u_name)
    if n_pass == n_tests
        fprintf('PASS')
    else
        fprintf('FAIL')
    end
    fprintf(': %d/%d tests pass (%g%%)\n\n', n_pass, n_tests, round(n_pass/n_tests*100,1))
    
    total_unittests = total_unittests + n_tests;
    total_unitpass = total_unitpass + n_pass;
end
unittest_time = toc(unittest_time);

%% ====== DOCTESTS =======
% Warn if using an older MATLAB than R2016b. Note version 9.1 = R2016b
if verLessThan('matlab','9.1')
    warning('dynamo_run_all_tests:matlab_version_too_old', ...
        'MATLAB changed its output formatting in R2016b. Many of these doctests will break on older versions')
end

fprintf('\n------ DOCTESTS -------\n')
docs_to_test = {    
                    %Sets
                    'setSingleItem'
                    'setList'
                    'setBasic'
                    'setCombinWithLimits'
                    %Random Processes
                    'rpDiscreteSample'
                    'rpLattice'
                    'rpTransMatrix'
                    %Utilities
                    'utilRandSetCurState'
                    'utilRandStatefromState'
                    'IntegerRangeFromReal'
                    'utilExtendRowVector'
                    %Examples
                    'doctest_MultiInv'
                    'MultiInv_demo'
                    'MultiInvSetupProblem'
                };
%initialize pass & test counts
total_doctests = 0;
total_docpass = 0;

doctest_time = tic;
for d_idx = 1:length(docs_to_test)
    d_name = docs_to_test{d_idx};
    
    fprintf('\nTesting %s\n', d_name)
    [~, n_pass, n_tests] = doctest(d_name);
    total_doctests = total_doctests + n_tests;
    total_docpass = total_docpass + n_pass;
end
doctest_time = toc(doctest_time);

%% Display summary
fprintf('\n========= dynamo testing summary ========\n')
%unit tests
fprintf(' UNIT TEST (%2d files): ', length(unittest_list))
if total_unitpass == total_unittests
    fprintf('PASS')
else
    fprintf('FAIL')
end
fprintf(', %3d/%3d suites pass (%g%%) in %g sec\n', total_unitpass, total_unittests, round(total_unitpass/total_unittests*100, 1), unittest_time)


%doctests
fprintf(' DOCTEST   (%2d files): ', length(docs_to_test))
if total_docpass == total_doctests
    fprintf('PASS')
else
    fprintf('FAIL')
end
fprintf(', %3d/%3d lines  pass (%g%%) in %g sec\n', total_docpass, total_doctests, round(total_docpass/total_doctests*100, 1), doctest_time)
end