function run_all_tests
% Run all tests for dynamo
%  Add additional tests, doctest or otherwise, as needed

docs_to_test = {    
                    %Sets
                    'setSingleItem'
                    %Random Processes
                    'rpDiscreteSample'
                    %Utilities
                    'utilRandSetCurState'
                    %Examples
                    'doctest_MultiInv'
                };

for d_idx = 1:length(docs_to_test)
    d_name = docs_to_test{d_idx};
    
    fprintf('Testing %s\n', d_name)
    doctest(d_name)
end
