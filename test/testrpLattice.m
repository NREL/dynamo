function tests = testrpLattice
%TESTrpLattice Test functions for rpLattice class
%
%  To use run the tests do the following:
%   1) Change the working directory to the folder that contains this file
%   2) enter the command: runtests at the command prompt (or in the
%   profiler)
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   8  2017-07-18 11:28  BryanP      BUGFIX: corrected inconsistant order in dlistnext to (t,s) 
%   7  2017-07-15 14:40  BryanP      Overhaulled tests for rpLattice v15 
%   6  2017-07-14 22:55  BryanP      Updated to use internal MATLAB unit test framework (post R2013) 
%   5  2012-04-17 08:55  BryanP      Reworked cdf for sim/sample 
%   4  2011-06-09 16:15  BryanP      Return to default RandStream, convert to get/setGlobalStream  
%   3  2011-01-04 22:00  BryanP      Complete, working version
%   2  2010-12-23 15:45  BryanP      Continued expansion & test for dlist*
%   1  2010-12-13 23:30  BryanP      Adapted from testrpList v3

    tests = functiontests(localfunctions);
end

%% Shared Setup, used by all tests
function setup(testCase)

    testCase.TestData.baseline_rand = RandStream.getGlobalStream;
    testCase.TestData.n_samples = 2e4;
    testCase.TestData.sample_tol = 0.02;    %Relative tolerance for mean of samples vs analytic from PMF

    testCase.TestData.start = 5;
    testCase.TestData.coef = [1/2  6/5 3/2 ]';  %Note: This is only partially geometric so some but not all states will recombine
    testCase.TestData.prob = [0.15 0.5 0.35]';
    testCase.TestData.tmax = 4;
    
    %Semi-hand calculated correct answer for lattice
    %Note: use l=start; then repeat l = unique(c * l'); rats(l)
    testCase.TestData.lattice = {5, ... 
                    [5/2 6   15/2]', ...
                    [5/4 3   15/4 36/5 9   45/4]', ...
                    [5/8 3/2 15/8 18/5 9/2 45/8 216/25 54/5 27/2 135/8]'};
    testCase.TestData.value_list = [5/8 5/4 3/2 15/8 5/2 3 18/5 15/4 9/2 5 45/8 6 36/5 15/2 216/25 9 54/5 45/4 27/2 135/8]';

	%Build up unconditional probability
    uncond_prob{1} = 1;
    uncond_prob{2} = testCase.TestData.prob;
    uncond_prob{3} = ...
        ... % [5/4  3    15/4 36/5 9    45/4]      %state(t+1)
              [0.15 0.5  0.35 0    0    0   ]' * uncond_prob{2}(1) ... %prob of state(t+1) given 5/2
            + [0    0.15 0    0.5  0.35 0   ]' * uncond_prob{2}(2) ... %prob of state(t+1) given 6
            + [0    0    0.15 0    0.5  0.35]' * uncond_prob{2}(3);    %prob of state(t+1) given 15/2
    %  normalize to 1
    uncond_prob{3} = uncond_prob{3}/sum(uncond_prob{3});
    uncond_prob{4} = ...
        ... % [5/8  3/2  15/8 18/5 9/2  45/8 216/25 54/5 27/2 135/8]      %state(t+1)
              [0.15 0.5  0.35 0    0    0    0      0    0    0    ]' * uncond_prob{3}(1) ... %prob of state(t+1) given 5/4
            + [0    0.15 0    0.5  0.35 0    0      0    0    0    ]' * uncond_prob{3}(2) ... %prob of state(t+1) given 3
            + [0    0    0.15 0    0.5  0.35 0      0    0    0    ]' * uncond_prob{3}(3) ... %prob of state(t+1) given 15/4
            + [0    0    0    0.15 0    0    0.5    0.35 0    0    ]' * uncond_prob{3}(4) ... %prob of state(t+1) given 36/5
            + [0    0    0    0    0.15 0    0      0.5  0.35 0    ]' * uncond_prob{3}(5) ... %prob of state(t+1) given 9
            + [0    0    0    0    0    0.15 0      0    0.5  0.35 ]' * uncond_prob{3}(6) ;   %prob of state(t+1) given 45/4
    %  normalize to 1
    uncond_prob{4} = uncond_prob{4}/sum(uncond_prob{4});
    %Store results
    testCase.TestData.uncond_prob = uncond_prob;
                
	%Actually build the demo object to use in other tests
    testCase.TestData.lattice_object = rpLattice(testCase.TestData.start, testCase.TestData.coef, testCase.TestData.prob, testCase.TestData.tmax);
end

%% Shared Teardown, used by all tests
function teardown(testCase)
    RandStream.setGlobalStream(testCase.TestData.baseline_rand);
end

%% Test constructor
function testConstructor(testCase)
    %Test basic constructor
    lattice_object = testCase.TestData.lattice_object;
    
    testCase.verifyEqual(lattice_object.Start, testCase.TestData.start)
    testCase.verifyEqual(lattice_object.Coef,  testCase.TestData.coef)
    testCase.verifyEqual(lattice_object.CondProb,  testCase.TestData.prob)
    testCase.verifyEqual(lattice_object.t, 1)
    
    %Ensure we convert row vectors to column vectors, and start with 0 for tmax:
    b = rpLattice(4, [2 3], ones(2,1)/2);
    testCase.verifyEqual(b.Coef,  [2 3]')
    testCase.verifyEqual(b.CondProb,  [0.5 0.5]')
    clear b
    
    b = rpLattice(4, [2 3]', ones(1,2)/2);
    testCase.verifyEqual(b.Coef,  [2 3]')
    testCase.verifyEqual(b.CondProb,  [0.5 0.5]')
    clear b

    
    %--Check for other constructor errors
    % coef & prob length mis-match
    bad = @() rpLattice(4, [2 3]', ones(3,1)/3);
    testCase.verifyError(bad, 'rpLattice:CoefProbMismatch');
    
    % non-scalar start value
    bad = @() rpLattice([4 2], [0.3 2 3]', ones(3,1)/3);
    testCase.verifyError(bad, 'rpLattice:NonScalarStart');

    % negative coeficients
    bad = @() rpLattice(4, [0.3 -2 3]', ones(3,1)/3);
    testCase.verifyError(bad, 'rpLattice:NegCoef');

    % sum(prob)>1
    bad = @() rpLattice(4, [0.3 2 3]', ones(3,1));
    testCase.verifyError(bad, 'rpLattice:ProbNotSumOne');
end

%% Test setparams method
function testSetParams(testCase)
    %Test basic function
    s = 9;
    c = [0.3 1 3]';
    p = ones(3,1)/3;
    t_max = 10;
    testCase.TestData.lattice_object.setparams(s, c, p, t_max);
    testCase.verifyEqual(testCase.TestData.lattice_object.Start, s)
    testCase.verifyEqual(testCase.TestData.lattice_object.Coef,  c)
    testCase.verifyEqual(testCase.TestData.lattice_object.CondProb,  p)
    
    %Ensure we convert row vectors to column vectors:
    testCase.TestData.lattice_object.setparams(4, [2 3], ones(2,1)/2);
    testCase.verifyEqual(testCase.TestData.lattice_object.Coef,  [2 3]')
    testCase.verifyEqual(testCase.TestData.lattice_object.CondProb,  [0.5 0.5]')
    
    testCase.TestData.lattice_object.setparams(4, [2 3]', ones(1,2)/2);
    testCase.verifyEqual(testCase.TestData.lattice_object.Coef,  [2 3]')
    testCase.verifyEqual(testCase.TestData.lattice_object.CondProb,  [0.5 0.5]')

    %--Check for error handling
    % coef & prob length mis-match
    bad = @() testCase.TestData.lattice_object.setparams(4, [2 3]', ones(3,1)/3);
    testCase.verifyError(bad, 'rpLattice:CoefProbMismatch');
    
    % non-scalar start value
    bad = @() testCase.TestData.lattice_object.setparams([4 2], [0.3 2 3]', ones(3,1)/3);
    testCase.verifyError(bad, 'rpLattice:NonScalarStart');

    % negative coeficients
    bad = @() testCase.TestData.lattice_object.setparams(4, [0.3 -2 3]', ones(3,1)/3);
    testCase.verifyError(bad, 'rpLattice:NegCoef');

    % sum(prob)>1
    bad = @() testCase.TestData.lattice_object.setparams(4, [0.3 2 3]', ones(3,1));
    testCase.verifyError(bad, 'rpLattice:ProbNotSumOne');
end

%% Test Setting individual parameters
function testIndividualParams(testCase)
    %New values (all same size so should work fine)
    s = testCase.TestData.start *3;
    c = (1:length(testCase.TestData.coef))';
    p = ones(size(testCase.TestData.prob))/length(testCase.TestData.prob);
    
    %Try these out
    testCase.TestData.lattice_object.Start = s;
    testCase.verifyEqual(testCase.TestData.lattice_object.Start, s)
    testCase.TestData.lattice_object.Coef = c;
    testCase.verifyEqual(testCase.TestData.lattice_object.Coef,  c)
    testCase.TestData.lattice_object.CondProb = p;
    testCase.verifyEqual(testCase.TestData.lattice_object.CondProb,  p)
        
    %--Check for error handling
    % coef & prob length mis-match
    try
    	testCase.TestData.lattice_object.Coef = [c; 42];
        error('test:junk', 'Previous Error not thrown')
    catch exception
        if not(strcmp(exception.identifier,'rpLattice:CoefProbMismatch'))
            throw(exception);
        end
    end
    
    %Set back to correct value to avoid problems
    testCase.TestData.lattice_object.Coef = c;
    try
        testCase.TestData.lattice_object.CondProb = [p; 0];
        error('test:junk', 'Previous Error not thrown')
    catch exception
        if not(strcmp(exception.identifier,'rpLattice:CoefProbMismatch'))
            throw(exception);
        end
    end
    
    testCase.TestData.lattice_object.CondProb = p;
    % non-scalar start value
    try
        testCase.TestData.lattice_object.Start = [4; 2];
        error('test:junk', 'Previous Error not thrown')
    catch exception
        if not(strcmp(exception.identifier,'rpLattice:NonScalarStart'))
            throw(exception);
        end
    end

    % negative coeficients
    c(1) = -c(1);
    try
        testCase.TestData.lattice_object.Coef = c;
        error('test:junk', 'Previous Error not thrown')
    catch exception
        if not(strcmp(exception.identifier, 'rpLattice:NegCoef'))
            throw(exception);
        end
    end

    % sum(prob)>1
    try
        testCase.TestData.lattice_object.CondProb = p * 2;
        error('test:junk', 'Previous Error not thrown')
    catch exception
        if not(strcmp(exception.identifier, 'rpLattice:ProbNotSumOne'))
            throw(exception);
        end
    end
end

%% Test as_array method
function testAsArray(testCase)
    % Test proper handling of t<1 
    bad = @() testCase.TestData.lattice_object.as_array(0);
    testCase.verifyError(bad, 'RandProcess:InvalidTime')
    
    % --Test proper handling of specified times for 1<=t<=Tmax+5
    % For all time periods
    for t = 1:testCase.TestData.tmax+5
        [v_list, p_list] = testCase.TestData.lattice_object.as_array(t);

        if t > testCase.TestData.tmax
            t_lookup = testCase.TestData.tmax;
        else
            t_lookup = t;
        end
        %Check outputs
        testCase.verifyEqual(length(v_list), length(p_list))
        testCase.verifyEqual(v_list, testCase.TestData.lattice{t_lookup})
        testCase.verifyEqual(p_list, testCase.TestData.uncond_prob{t_lookup})
        
        %Verify that object time has not advanced
        testCase.verifyEqual(testCase.TestData.lattice_object.t, 1)
    end
    
    % Test reduced number of outputs
end

%% Test sample method
%
% Note: the conditional sampling portion (nargin=4) extensively tested in
% step()
function testSample(testCase)
    rng(0) %Reset random number generator to ensure consistant results for samples
    
    % Test proper handling of t<1 
    bad = @() testCase.TestData.lattice_object.sample([], 0);
    testCase.verifyError(bad, 'RandProcess:InvalidTime')
    
    %Verify samples over time
    for t = 1:testCase.TestData.tmax+2
        s_list = testCase.TestData.lattice_object.sample(testCase.TestData.n_samples, t);
        
        if t > testCase.TestData.tmax
            t_lookup = testCase.TestData.tmax;
        else
            t_lookup = t;
        end
        % Check correct # samples
        testCase.verifyEqual(size(s_list,1), testCase.TestData.n_samples);
        % Check that all values are valid for this timeperiod
        testCase.verifyTrue(all(ismember(s_list, testCase.TestData.lattice{t_lookup})));
        % Check sample mean is within relative tolerance of analytic PMF
        testCase.verifyEqual(mean(s_list), testCase.TestData.lattice{t_lookup}' * testCase.TestData.uncond_prob{t_lookup}, 'RelTol', testCase.TestData.sample_tol)
        
        %Verify that object time has not advanced
        testCase.verifyEqual(testCase.TestData.lattice_object.t, 1)
    end
        
    % Test reduced number of outputs
end

%% Test Dlistprev method
function testDlistPrev(testCase)
    % Test proper handling of t<0
    bad = @() testCase.TestData.lattice_object.dlistprev(-1, testCase.TestData.lattice{1});
    testCase.verifyError(bad, 'RandProcess:InvalidTime')
    
    % Test proper handling of t=0
    bad = @() testCase.TestData.lattice_object.dlistprev(0, testCase.TestData.lattice{1});
    testCase.verifyError(bad, 'RandProcess:InvalidTime')

    % Test proper handling of t=1 (because can't go backwards from t=1
    bad = @() testCase.TestData.lattice_object.dlistprev(1, testCase.TestData.lattice{1});
    testCase.verifyError(bad, 'RandProcess:InvalidTime')

    % --Test proper handling of previous states for t = 2<=t<=Tmax
    % For all time periods
    for t = 2:testCase.TestData.tmax
        %Valid States & times
        CheckPrevStates(testCase, t)
    end

    % Test for constant output for t>Tmax
    t_max = testCase.TestData.tmax;
    for s_idx = 1:length(testCase.TestData.lattice{t_max})
        s=testCase.TestData.lattice{t_max}(s_idx);

        for t = t_max+1:t_max+10
            [prev_v, p] = testCase.TestData.lattice_object.dlistprev(t, s);
            testCase.verifyEqual(prev_v, s);
            testCase.verifyEqual(p, 1);
        end
    end
    
    
    %Test a totally invalid state
    bad = @() testCase.TestData.lattice_object.dlistprev(2, 98);
    testCase.verifyError(bad, 'RandProcess:InvalidState')
        
    %Test using testCase.TestData.lattice_object valid states but at the wrong time
    for t = 2:testCase.TestData.tmax
        bad = @() testCase.TestData.lattice_object.dlistprev(t, testCase.TestData.lattice{t-1}(1));
        testCase.verifyError(bad, 'RandProcess:InvalidState')
    end
    
    
    % Test use of current (unspecified) time
    % t=2;
    testCase.TestData.lattice_object.step();
    CheckPrevStates(testCase)
    
    % Test reduced number of outputs
    % t = 3;
    testCase.TestData.lattice_object.step();
    % value only
    CheckPrevStates(testCase, t, 1)
    % value & probablity only
    CheckPrevStates(testCase, t, 2)
end

%% Test Dlistnext method
function testDlistNext(testCase)
    % Test proper handling of t<-1 
    bad = @() testCase.TestData.lattice_object.dlistnext(-2, testCase.TestData.lattice{1});
    testCase.verifyError(bad, 'RandProcess:InvalidTime')
    
    % --Test proper handling of -1<=t<=Tmax-1
    % For all time periods
    for t = 1:testCase.TestData.tmax-1
        % For all possible states in this time period
        for s_idx = 1:length(testCase.TestData.lattice{t})
            s=testCase.TestData.lattice{t}(s_idx);
            [next_v, p] = testCase.TestData.lattice_object.dlistnext(t,s);
            
            testCase.verifyEqual(next_v, s .* testCase.TestData.coef, 'AbsTol', 1e-6)
            testCase.verifyEqual(p, testCase.TestData.prob)
        end
    end
    
    % Test for constant output for t>Tmax
    t_max = testCase.TestData.tmax;
    for s_idx = 1:length(testCase.TestData.lattice{t_max})
        s=testCase.TestData.lattice{t_max}(s_idx);

        for t = t_max+1:t_max+10
            [next_v, p] = testCase.TestData.lattice_object.dlistnext(t,s);
            testCase.verifyEqual(next_v, s);
            testCase.verifyEqual(p, 1);
        end
    end
    
    %Test a totally invalid state
    bad = @() testCase.TestData.lattice_object.dlistnext(1, 98);
    testCase.verifyError(bad, 'RandProcess:InvalidState')
        
    %Test using a valid states but at the wrong time
    for t = 1:testCase.TestData.tmax-1
        bad = @() testCase.TestData.lattice_object.dlistnext(t, testCase.TestData.lattice{t+1}(1));
        testCase.verifyError(bad, 'RandProcess:InvalidState')
    end
    
    % Test use of current (defaults to t=1) time
    t=1;
    for s_idx = 1:length(testCase.TestData.lattice{t})
        s=testCase.TestData.lattice{t}(s_idx);
        [next_v, p] = testCase.TestData.lattice_object.dlistnext([], s);

        testCase.verifyEqual(next_v, RoundTo(s .* testCase.TestData.coef, testCase.TestData.lattice_object.Tol))
        testCase.verifyEqual(p, testCase.TestData.prob)
    end

    % Test reduced number of outputs
    t = 2;
    testCase.TestData.lattice_object.step();
    % value only
    for s_idx = 1:length(testCase.TestData.lattice{t})
        s=testCase.TestData.lattice{t}(s_idx);
        [next_v] = testCase.TestData.lattice_object.dlistnext(t,s);

        testCase.verifyEqual(next_v, RoundTo(s .* testCase.TestData.coef, testCase.TestData.lattice_object.Tol))
    end   
end

%% Test Range method
function testRange(testCase)
    % Test extraction of all possible states
    v = testCase.TestData.lattice_object.range('all');
    testCase.verifyEqual(v, testCase.TestData.value_list([1, end]))
    
    % Test proper handling of t<1
    bad = @() testCase.TestData.lattice_object.range(-0.1);
    testCase.verifyError(bad, 'RandProcess:InvalidTime');

    bad = @() testCase.TestData.lattice_object.range(0);
    testCase.verifyError(bad, 'RandProcess:InvalidTime');

    % Test proper handling of 1<t<Tmax
    for t = 1:testCase.TestData.tmax
        v = testCase.TestData.lattice_object.range(t);
        testCase.verifyEqual(v, [testCase.TestData.lattice{t}(1); testCase.TestData.lattice{t}(end)]);
    end
    
    % Test proper handling of t>Tmax
    t_max = testCase.TestData.tmax;
    for t = t_max+1:t_max+10
        v = testCase.TestData.lattice_object.range(t);
        testCase.verifyEqual(v, [testCase.TestData.lattice{t_max}(1); testCase.TestData.lattice{t_max}(end)]);
    end
    
    % Test use of current (unspecified) time
    t=2;
    
    testCase.TestData.lattice_object.step();
    
    v = testCase.TestData.lattice_object.range();
    testCase.verifyEqual(v, [testCase.TestData.lattice{t}(1); testCase.TestData.lattice{t}(end)]);
        
    % Test proper handling of non-integer times
    for t = [1.2 4.1 3.1 pi]
        v = testCase.TestData.lattice_object.range(t);
        t = floor(t); %#ok<FXSET>
        testCase.verifyEqual(v, [testCase.TestData.lattice{t}(1); testCase.TestData.lattice{t}(end)]);
    end
    
    % Test reduced number of outputs
    t = 3;
    % value only
    testCase.TestData.lattice_object.step();
    v = testCase.TestData.lattice_object.range(t);
    testCase.verifyEqual(v, [testCase.TestData.lattice{t}(1); testCase.TestData.lattice{t}(end)]);
    
end

%% Test Step method
function testStep(testCase)
    %-- Test proper action when steping through the range of t
    % Note: the internal t starts at 1

    t_max = testCase.TestData.tmax + 10;
    % Find "right" random based answers, including t>Tmax
    v_test = SimHelper(testCase, 1, t_max, testCase.TestData.start);  

    for t = 2:t_max %Start at 2 b/c comparing results after step()
        [v, t_out] = testCase.TestData.lattice_object.step();
        testCase.verifyEqual(v, v_test(t));
        testCase.verifyEqual(testCase.TestData.lattice_object.t, t)
        testCase.verifyEqual(t, t_out)
        testCase.verifyEqual(testCase.TestData.lattice_object.cur_state(), v_test(t))
    end

    %-- Test bigger steps
    t = 1;

    for st = 1:3
        % First reset since that may change random samples
        testCase.TestData.lattice_object.reset();

        % Find "right" random based answers
        v_test = SimHelper(testCase, t, testCase.TestData.tmax+1, testCase.TestData.start);  

        idx = st+1;
        [v, t_out] = testCase.TestData.lattice_object.step(st);
        testCase.verifyEqual(v, v_test(idx));
        testCase.verifyEqual(testCase.TestData.lattice_object.t, t+st)
        testCase.verifyEqual(t+st, t_out)
    end

    % Test zero stepsize (should have no change)
    t = testCase.TestData.lattice_object.t;
    st = 0;
    [v, t_out] = testCase.TestData.lattice_object.step(st);
    testCase.verifyEqual(v, v_test(idx));
    testCase.verifyEqual(testCase.TestData.lattice_object.t, t+st)
    testCase.verifyEqual(t+st, t_out)
    
    %test non-integer steps
    t = 1;
    for st = [1.2 0.1 2 pi-1]
        % First reset since that may change random samples
        testCase.TestData.lattice_object.reset();

        % Find "right" random based answers
        v_test = SimHelper(testCase, t, testCase.TestData.tmax+1, testCase.TestData.start);  

        idx = floor(st)+1;
        [v, t_out] = testCase.TestData.lattice_object.step(st);
        testCase.verifyEqual(v, v_test(idx));
        testCase.verifyEqual(testCase.TestData.lattice_object.t, floor(t+st))
        testCase.verifyEqual(floor(t+st), t_out)    
    end

    % Test reduced number of outputs & starting with testCase.TestData.lattice_object non-zero start time
    t = 1;
    % Find "right" random based answers
    v_test = SimHelper(testCase, t, testCase.TestData.tmax, testCase.TestData.start);  

    testCase.TestData.lattice_object.reset();
    % value only
    v = testCase.TestData.lattice_object.step();
    testCase.verifyEqual(v, v_test(2));
    % value & state only
    v = testCase.TestData.lattice_object.step();
    testCase.verifyEqual(v, v_test(3));
    
end

%% Test CurState method
function testCurState(testCase)

    v_test = SimHelper(testCase, 1, testCase.TestData.tmax+5, testCase.TestData.start);  

    % Test proper handling of 1<t<Tmax+
    for t = 1:(length(testCase.TestData.tmax+4))
        testCase.TestData.lattice_object.step();
        
        v = testCase.TestData.lattice_object.cur_state();
        testCase.verifyEqual(v, v_test(t+1));
    end 
end

%% Test Reset method
function testReset(testCase)
    testCase.TestData.lattice_object.reset();
    
    testCase.verifyEqual(testCase.TestData.lattice_object.t, 1);
    testCase.verifyEqual(testCase.TestData.lattice_object.cur_state(), testCase.TestData.start);
end





%% ---------------------------- Helper Functions -------------------------------
function CheckPrevStates(testCase, t, num_out)
    if nargin < 2
        t = testCase.TestData.lattice_object.t;
        specify_t = false;
    else
        specify_t = true;
    end
    
    if nargin < 3
        num_out = 2;
    end
        
    
    % For all possible states in this time period
    for s_idx = 1:length(testCase.TestData.lattice{t})
        s=testCase.TestData.lattice{t}(s_idx);
        if specify_t
            if num_out == 2
                [prev_v, p] = testCase.TestData.lattice_object.dlistprev(t,s);
            else
                prev_v = testCase.TestData.lattice_object.dlistprev(t,s);
            end
        else
            if num_out == 2
                [prev_v, p] = testCase.TestData.lattice_object.dlistprev([],s);
            else
                prev_v = testCase.TestData.lattice_object.dlistprev([],s);
            end
        end

        % make sure that all possible prior states are listed with the
        % correct probabilities by looping over all transition
        % coeficients
        vals_found = 0;
        for c_idx = 1:length(testCase.TestData.coef)
            possible_v = RoundTo(s/testCase.TestData.coef(c_idx), 0.0001);
            [found_v, ~] = intersect(possible_v, testCase.TestData.lattice{t}); 
            if isempty(found_v)
                continue
            elseif length(found_v) > 1
                %If we get here, our lattice is wrong
                error('testrpLattice:toomanyprev', 'Multiple previous states found for one coefficient')
            else
                vals_found = vals_found+1;
                testCase.verifyTrue(ismember(found_v, prev_v))
%TODO: Check for proper probabilities
            end
            testCase.verifyTrue(vals_found>0)
        end
        if num_out == 2
            testCase.verifyEqual(sum(p), 1, 'AbsTol', 1e-8)
        end
    end
end

function v_test = SimHelper(testCase, t_min, t_max, cur_val, seed)
    if nargin < 5 || isempty(seed)
        c = clock;          %array of time value
        seed = c(end)*1000; %Use current seconds of time in ms
    end
    
    % reset the random stream using a seed-based rand sequence
    % Note: we assume the rands are called in simulation order and that all
    % rand values are used. This may have to change with alternate
    % implementations
    rs=RandStream('mt19937ar', 'Seed', seed);
    RandStream.setGlobalStream(rs);
    
    v_test = NaN(t_max-t_min+1, 1);
    v_test(1) = cur_val;
    
    for t_step = 2: (t_max-t_min+1)
        if t_step > testCase.TestData.tmax
            v_test(t_step) = v_test(t_step-1);
        else
            next = find(rand <= cumsum(testCase.TestData.prob), 1, 'first');
            v_test(t_step) = RoundTo(v_test(t_step-1)*testCase.TestData.coef(next), 0.0001);
        end
    end

    % reset the random stream using a seed-based rand sequence
    rs=RandStream('mt19937ar', 'Seed', seed);
    RandStream.setGlobalStream(rs);

end
