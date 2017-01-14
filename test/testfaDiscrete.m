function test_suite = testfaDiscrete %#ok<MCUOA>
%TESTfaDiscrete Test functions for faDiscrete class
%
%  To use run the tests do the following:
%   1) Change the working directory to the folder that contains this file
%   2) enter the command: runtests at the command prompt (or in the
%   profiler)
%
%   Note: this code relies on MATLAB xUnit from the MATLAB file exchange.
%   You can get it at: 
%    http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework
%   Be sure to add the xunit folder to your MATLAB path
%
% TODO: EXPAND!! the current version provides only very basic checks
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-03-24 15:30  BryanP      Adapted from testrpLattice v3

    initTestSuite;
end

%% Shared Setup, used by all tests
function a = setup %#ok<*DEFNU>
    global TfaD_max_per_dim
    
    TfaD_max_per_dim = [1 2 3];

	%Actually build the demo object to use in other tests
    a = faDiscrete(TfaD_max_per_dim);
end

%% Shared Teardown, used by all tests
function teardown(a)
    global TfaD_max_per_dim %#ok<NUSED>

    a.delete;
    clearvars a TfaD*
end

%% Test constructor
function testConstructor(a)
    global TfaD_max_per_dim

    %Basic loop to test different constructors
    for count = 1:2
        %Test basic constructor
        assertEqual(a.N_dim, length(TfaD_max_per_dim))
        assertEqual(a.Range, vertcat(ones(size(TfaD_max_per_dim)), TfaD_max_per_dim))

        assertEqual(a.Needs_jacobian, false)
        
        CheckAllStates(a, zeros(TfaD_max_per_dim))
    
        %Ensure we convert column vectors to row vectors:
        %
        % Note: this uses the new format for the second pass
        a = faDiscrete(TfaD_max_per_dim);
    end
end

%% Test update & approx methods
%
% Note: assumes we are using 1/n stepsize (default)
function testUpdateApprox(a)

    % Simple one value
    [v, ~, step] = a.update([1 2 1], 1);
    assertEqual(v, 1);
    assertEqual(step, 1);
    
    [v, ~, step] = a.approx([1 2 1]);
    assertEqual(v, 1);
    assertEqual(step, 1);

    % Update more than one value
    [v, ~, step] = a.update([1 1 1; 1 1 2; 1 2 1], [2;5;5]);
    assertEqual(v, [2;5;3]);
    assertEqual(step, [1;1;0.5]);

end

%% Test finding Max
function testMax(a)

    % Run a few updates
    a.update([1 2 2], 10);
    a.update([1 1 3; 1 1 2; 1 2 2], [2;4;5]);
    a.update([1 2 3; 1 1 1], [6;5]);
    
    % Note at this point our approximation looks like:
    % (all with first index = 1):
    % [5   4 2
    %  0 7.5 6]

    %find the unconstrained maximum
    [max_val, state] = a.max();
    assertEqual(max_val, 7.5)   %(10 + 5)/2
    assertEqual(state, [1 2 2])
    
    %max_limit constrained maximum
    [max_val, state] = a.max([1 1 3]);
    assertEqual(max_val, 5)
    assertEqual(state, [1 1 1])
    
    %min_limit (only) constrained maximum
    [max_val, state] = a.max([], [1 1 3]);
    assertEqual(max_val, 6)
    assertEqual(state, [1 2 3])

    %min & max limit constrained maximum
    [max_val, state] = a.max([1 1 3], [1 1 2]);
    assertEqual(max_val, 4)
    assertEqual(state, [1 1 2])

end


%% Test Max Scaled Sum
function testMaxScaledSum(a)

    % Run a few updates
    a.update([1 2 2], 10);
    a.update([1 1 3; 1 1 2; 1 2 2], [2;4;5]);
    a.update([1 2 3; 1 1 1], [6;5]);
    
    % Note at this point our approximation looks like:
    % (all with first index = 1):
    % [5   4  2
    %  0 7.5  6]

    %-- First consider adding to ourselves with a 0.5 & 2 scaling
    % Since the scaling is positive the -100 in b will never be choosen and
    % our results will simply be scaled versions of the max() results

    b = a.copy();
    b.update([1 1 3], -100);
    % Note at this point b looks like:
    % (all with first index = 1):
    % [5   4  -49
    %  0 7.5    6]
    
    
    %find the unconstrained maximum
    [max_val, state] = a.maxScaledSum(b, [0.5 2]);
    scale_factor = 2.5;
    assertEqual(max_val, 7.5*scale_factor)   %(10 + 5)/2
    assertEqual(state, [1 2 2])
    
    %max_limit constrained maximum
    [max_val, state] = a.maxScaledSum(b, [0.5 2], [1 1 3]);
    assertEqual(max_val, 5*scale_factor)
    assertEqual(state, [1 1 1])
    
    %min_limit (only) constrained maximum
    [max_val, state] = a.maxScaledSum(b, [0.5 2], [], [1 1 3]);
    assertEqual(max_val, 6*scale_factor)
    assertEqual(state, [1 2 3])

    %min & max limit constrained maximum
    [max_val, state] = a.maxScaledSum(b, [0.5 2], [1 1 3], [1 1 2]);
    assertEqual(max_val, 4*scale_factor)
    assertEqual(state, [1 1 2])

    
    %-- Now make that negative the max with [1 -1] scaling
    [max_val, state] = a.maxScaledSum(b, [1 -1]);
    assertEqual(max_val, 49+2)
    assertEqual(state, [1 1 3])
    
    
    %make sure we don't have any lingering objects
    clear b
end


%% ---------------------------- Helper Functions -------------------------------

% Recursively compare the current function approximation with provided data
function CheckAllStates(a, data, dim_max, idx_list)
    global TfaD_max_per_dim

    % If called with no arguements, setup additional parameters required for
    % recursive use
    if nargin < 3
        dim_max = TfaD_max_per_dim;
        idx_list = [];
    end
    
    % Loop over the range of the first dimension
    for idx = 1:dim_max(1)
        
        %build up index list for this iteration
        this_idx_list = [idx_list idx];
        
        if isscalar(dim_max)
            %recursive base case: 
            % check that our value matches the approximation
            assertEqual(a.approx(this_idx_list), ...
                        data(mat2ind(size(data), this_idx_list)) )
        else
            %Recursive call:
            % Use our current iteration for the next index and adjust dim_max
            CheckAllStates(a, data, dim_max(2:end), this_idx_list)
        end
    end
    
end

