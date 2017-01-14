function test_suite = testrpLattice %#ok<MCUOA>
%TESTrpLattice Test functions for rpLattice class
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
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-12-13 23:30  BryanP      Adapted from testrpList v3
%   2  2010-12-23 15:45  BryanP      Continued expansion & test for dlist*
%   3  2011-01-04 22:00  BryanP      Complete, working version
%   4  2011-06-09 16:15  BryanP      Return to default RandStream, convert to get/setGlobalStream  
%   5  2012-04-17 08:55  BryanP      Reworked cdf for sim/sample 


    initTestSuite;
end

%% Shared Setup, used by all tests
function a = setup %#ok<*DEFNU>
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice 
    global TLtR_value_list TLtR_num_lattice
    global TLtR_baseline_rand
    TLtR_baseline_rand = RandStream.getGlobalStream;

    TLtR_start = 5;
    TLtR_coef = [0.5  1.2 1.5 ]';
    TLtR_prob = [0.15 0.5 0.35]';
    TLtR_tmax = 3;
    
    %Semi-hand calculated correct answer for lattice
    %Note: use l=start; then repeat l = unique(c * l'); rats(l)
    TLtR_lattice = {5, ... 
                    [5/2 6   15/2]', ...
                    [5/4 3   15/4 36/5 9   45/4]', ...
                    [5/8 3/2 15/8 18/5 9/2 45/8 216/25 54/5 27/2 135/8]'};
    % State numbers     1   2   3    4   5  6   7    8   9  10  11 12  13   14   15   16  17   18   19   20   
    TLtR_value_list = [5/8 5/4 3/2 15/8 5/2 3 18/5 15/4 9/2 5 45/8 6 36/5 15/2 216/25 9 54/5 45/4 27/2 135/8]';
    TLtR_num_lattice = { 10, ... 
                        [ 5 12 14]', ...
                        [ 2  6  8 13 16 18]', ...
                        [ 1  3  4  7  9 11 15 17 19 20]'};
                    
	%Actually build the demo object to use in other tests
    a = rpLattice(TLtR_start, TLtR_coef, TLtR_prob, TLtR_tmax);
end

%% Shared Teardown, used by all tests
function teardown(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice %#ok<NUSED>

    global TLtR_baseline_rand
    RandStream.setGlobalStream(TLtR_baseline_rand);

    a.delete;
    clearvars a TLtR_*
end

%% Test constructor
function testConstructor(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice %#ok<NUSED>

    %Test basic constructor
    assertEqual(a.Start, TLtR_start)
    assertEqual(a.Coef,  TLtR_coef)
    assertEqual(a.Prob,  TLtR_prob)
    assertEqual(a.t, 0)
    assertEqual(a.Tmax, TLtR_tmax)
    
    %Ensure we convert row vectors to column vectors, and start with 0 for tmax:
    b = rpLattice(4, [2 3], ones(2,1)/2);
    assertEqual(b.Coef,  [2 3]')
    assertEqual(b.Prob,  [0.5 0.5]')
    assertEqual(b.Tmax, 0)
    clear b
    
    b = rpLattice(4, [2 3]', ones(1,2)/2);
    assertEqual(b.Coef,  [2 3]')
    assertEqual(b.Prob,  [0.5 0.5]')
    clear b

    
    %--Check for other constructor errors
    % coef & prob length mis-match
    bad = @() rpLattice(4, [2 3]', ones(3,1)/3);
    assertExceptionThrown(bad, 'rpLattice:CoefProbMismatch');
    
    % non-scalar start value
    bad = @() rpLattice([4 2], [0.3 2 3]', ones(3,1)/3);
    assertExceptionThrown(bad, 'rpLattice:NonScalarStart');

    % negative coeficients
    bad = @() rpLattice(4, [0.3 -2 3]', ones(3,1)/3);
    assertExceptionThrown(bad, 'rpLattice:NegCoef');

    % sum(prob)>1
    bad = @() rpLattice(4, [0.3 2 3]', ones(3,1));
    assertExceptionThrown(bad, 'rpLattice:ProbNotSumOne');
end

%% Test setparams method
function testSetParams(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice %#ok<NUSED>

    %Test basic function
    s = 9;
    c = [0.3 1 3]';
    p = ones(3,1)/3;
    t_max = 10;
    a.setparams(s, c, p, t_max);
    assertEqual(a.Start, s)
    assertEqual(a.Coef,  c)
    assertEqual(a.Prob,  p)
    assertEqual(a.Tmax, 10)
    
    %Ensure we convert row vectors to column vectors:
    a.setparams(4, [2 3], ones(2,1)/2);
    assertEqual(a.Coef,  [2 3]')
    assertEqual(a.Prob,  [0.5 0.5]')
    
    a.setparams(4, [2 3]', ones(1,2)/2);
    assertEqual(a.Coef,  [2 3]')
    assertEqual(a.Prob,  [0.5 0.5]')

    %--Check for error handling
    % coef & prob length mis-match
    bad = @() a.setparams(4, [2 3]', ones(3,1)/3);
    assertExceptionThrown(bad, 'rpLattice:CoefProbMismatch');
    
    % non-scalar start value
    bad = @() a.setparams([4 2], [0.3 2 3]', ones(3,1)/3);
    assertExceptionThrown(bad, 'rpLattice:NonScalarStart');

    % negative coeficients
    bad = @() a.setparams(4, [0.3 -2 3]', ones(3,1)/3);
    assertExceptionThrown(bad, 'rpLattice:NegCoef');

    % sum(prob)>1
    bad = @() a.setparams(4, [0.3 2 3]', ones(3,1));
    assertExceptionThrown(bad, 'rpLattice:ProbNotSumOne');
end

%% Test Setting individual parameters
function testIndividualParams(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice %#ok<NUSED>

    %New values (all same size so should work fine)
    s = TLtR_start *3;
    c = (1:length(TLtR_coef))';
    p = ones(size(TLtR_prob))/length(TLtR_prob);
    t = 8;
    
    %Try these out
    a.Start = s;
    assertEqual(a.Start, s)
    a.Coef = c;
    assertEqual(a.Coef,  c)
    a.Prob = p;
    assertEqual(a.Prob,  p)
    a.Tmax = t;
    assertEqual(a.Tmax,  t)
        
    %--Check for error handling
    % coef & prob length mis-match
    try
    	a.Coef = [c; 42];
        error('test:junk', 'Previous Error not thrown')
    catch exception
        if not(strcmp(exception.identifier,'rpLattice:CoefProbMismatch'))
            throw(exception);
        end
    end
    
    %Set back to correct value to avoid problems
    a.Coef = c;
    try
        a.Prob = [p; 0];
        error('test:junk', 'Previous Error not thrown')
    catch exception
        if not(strcmp(exception.identifier,'rpLattice:CoefProbMismatch'))
            throw(exception);
        end
    end
    
    a.Prob = p;
    % non-scalar start value
    try
        a.Start = [4; 2];
        error('test:junk', 'Previous Error not thrown')
    catch exception
        if not(strcmp(exception.identifier,'rpLattice:NonScalarStart'))
            throw(exception);
        end
    end

    % negative coeficients
    c(1) = -c(1);
    try
        a.Coef = c;
        error('test:junk', 'Previous Error not thrown')
    catch exception
        if not(strcmp(exception.identifier, 'rpLattice:NegCoef'))
            throw(exception);
        end
    end

    % sum(prob)>1
    try
        a.Prob = p * 2;
        error('test:junk', 'Previous Error not thrown')
    catch exception
        if not(strcmp(exception.identifier, 'rpLattice:ProbNotSumOne'))
            throw(exception);
        end
    end
    
    % t > Tmax
    a.t = t + 1;
    assertEqual(a.t,  t+1)
   
    
    % t < 0
    try
        a.t = -1;
        error('test:junk', 'Previous Error not thrown')
    catch exception
        if not(strcmp(exception.identifier, 'RandProcess:InvalidTime'))
            throw(exception);
        end
    end
end

%% Test Dlist method
function testDlist(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice

    % Test extraction of all possible states
    [v, n] = a.dlist('all');
    assertEqual(v, TLtR_value_list)
    assertEqual(n, 1:length(TLtR_value_list))
    
    % Test proper handling of t<0
    bad = @() a.dlistprev(-1);
    assertExceptionThrown(bad, 'RandProcess:InvalidTime')
        
    % Test proper handling of 0<=t<=Tmax
    for t = 0:a.Tmax
        [v, n] = a.dlist(t);
        assertEqual(v, TLtR_lattice{t+1});
        assertEqual(n, TLtR_num_lattice{t+1});
    end
    
    % Test for constant output for t>Tmax
    t_max = a.Tmax;
    for t = t_max:t_max+10
        [v, n] = a.dlist(t);
        assertEqual(v, TLtR_lattice{t_max+1});
        assertEqual(n, TLtR_num_lattice{t_max+1});
    end
    
    % Test use of current (unspecified) time
    t=2;
    
    a.t = t;
    [v, n] = a.dlist();
    assertEqual(v, TLtR_lattice{t+1});
    assertEqual(n, TLtR_num_lattice{t+1});

    % Test reduced number of outputs
    t = 3;
    a.t = t;
    % value only
    v = a.dlist();
    assertEqual(v, TLtR_lattice{t+1});
end

%% Test Dlistprev method
function testDlistPrev(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice 

    % Test proper handling of t<0
    bad = @() a.dlistprev(TLtR_num_lattice{1},-1);
    assertExceptionThrown(bad, 'RandProcess:InvalidTime')
    
    % Test proper handling of t=0
    bad = @() a.dlistprev(TLtR_num_lattice{1},0);
    assertExceptionThrown(bad, 'RandProcess:InvalidTime')

    % --Test proper handling of 1<=t<=Tmax
    % For all time periods
    for t = 1:a.Tmax
        %Valid States & times
        CheckPrevStates(a, t)
    end

    % Test for constant output for t>Tmax
    t_max = a.Tmax;
    for s_idx = 1:length(TLtR_num_lattice{t_max+1})
        s=TLtR_num_lattice{t_max+1}(s_idx);
        this_v = TLtR_value_list(s);
        this_n = TLtR_num_lattice{t_max+1}(s_idx);

        for t = t_max+1:t_max+10
            [prev_v, prev_n, p] = a.dlistprev(s, t);
            assertEqual(prev_v, this_v);
            assertEqual(prev_n, this_n);
            assertEqual(p, 1);
        end
    end
    
    
    %Test a totally invalid state
    bad = @() a.dlistprev(98,1);
    assertExceptionThrown(bad, 'RandProcess:InvalidState')
        
    %Test using a valid states but at the wrong time
    for t = 1:a.Tmax
        bad = @() a.dlistprev(TLtR_num_lattice{t}(1),t);
        assertExceptionThrown(bad, 'RandProcess:InvalidState')
    end
    
    
    % Test use of current (unspecified) time
    t=2;
    
    a.t = t;
    CheckPrevStates(a)
    
    % Test reduced number of outputs
    t = 3;
    a.t = t;
    % value only
    CheckPrevStates(a, t, 1)
    % value & state only
    CheckPrevStates(a, t, 2)
end

%% Test Dlistnext method
function testDlistNext(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice

    % Test proper handling of t<-1 
    bad = @() a.dlistprev(TLtR_num_lattice{1},-2);
    assertExceptionThrown(bad, 'RandProcess:InvalidTime')
    
    % --Test proper handling of -1<=t<=Tmax-1
    % For all time periods
    for t = 0:a.Tmax-1
        % For all possible states in this time period
        for s_idx = 1:length(TLtR_num_lattice{t+1})
            s=TLtR_num_lattice{t+1}(s_idx);
            this_v = TLtR_value_list(s);
            [next_v, n, p] = a.dlistnext(s,t);
            
            assertEqual(next_v, RoundTo(this_v .* TLtR_coef, a.Tol))
            assertEqual(next_v, TLtR_value_list(n));
            assertEqual(p, TLtR_prob)
        end
    end
    
    % Test for constant output for t>Tmax
    t_max = a.Tmax;
    for s_idx = 1:length(TLtR_num_lattice{t_max+1})
        s=TLtR_num_lattice{t_max+1}(s_idx);
        this_v = TLtR_value_list(s);
        this_n = TLtR_num_lattice{t_max+1}(s_idx);

        for t = t_max:t_max+10
            [next_v, next_n, p] = a.dlistnext(s, t);
            assertEqual(next_v, this_v);
            assertEqual(next_n, this_n);
            assertEqual(p, 1);
        end
    end
    
    %Test a totally invalid state
    bad = @() a.dlistnext(98,1);
    assertExceptionThrown(bad, 'RandProcess:InvalidState')
        
    %Test using a valid states but at the wrong time
    for t = 0:a.Tmax-1
        bad = @() a.dlistnext(TLtR_num_lattice{t+2}(1),t);
        assertExceptionThrown(bad, 'RandProcess:InvalidState')
    end
    
    % Test use of current (unspecified) time
    t=1;
    
    a.t = t;
    for s_idx = 1:length(TLtR_num_lattice{t+1})
        s=TLtR_num_lattice{t+1}(s_idx);
        this_v = TLtR_value_list(s);
        [next_v, n, p] = a.dlistnext(s);

        assertEqual(next_v, RoundTo(this_v .* TLtR_coef, a.Tol))
        assertEqual(next_v, TLtR_value_list(n));
        assertEqual(p, TLtR_prob)
    end

    % Test reduced number of outputs
    t = 2;
    a.t = t;
    % value only
    for s_idx = 1:length(TLtR_num_lattice{t+1})
        s=TLtR_num_lattice{t+1}(s_idx);
        this_v = TLtR_value_list(s);
        [next_v] = a.dlistnext(s,t);

        assertEqual(next_v, RoundTo(this_v .* TLtR_coef, a.Tol))
    end
    
    % value & state only
    for s_idx = 1:length(TLtR_num_lattice{t+1})
        s=TLtR_num_lattice{t+1}(s_idx);
        this_v = TLtR_value_list(s);
        [next_v, n] = a.dlistnext(s,t);

        assertEqual(next_v, RoundTo(this_v .* TLtR_coef, a.Tol))
        assertEqual(next_v, TLtR_value_list(n));
    end
    
end


%% Test Dnum2val method
function testDnum2val(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice %#ok<NUSED>
    
    %test call with single, valid state number
    for idx = 1:length(TLtR_value_list)
        v = a.dnum2val(idx);
        assertEqual(v, TLtR_value_list(idx))
    end
    
    %test call with vector of state numbers
    idx = [1 3 1 5]';
    v = a.dnum2val(idx);
    assertEqual(v, TLtR_value_list(idx))
    
    %check that error thrown when passed a list including an invalid state num
    bad = @() a.dnum2val([2 10 52]);
    assertExceptionThrown(bad, 'rpLattice:InvalidStateNum');
end

%% Test Dval2num method
function testDval2num(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice %#ok<NUSED>
    
    %test call with single, valid values
    for idx = 1:length(TLtR_value_list)
        v = TLtR_value_list(idx);
        n = a.dval2num(v);
        assertEqual(n, idx)
    end
    
    %test call with vector of values
    idx = [1 4 1 5]';
    v = TLtR_value_list(idx);
    n = a.dval2num(v);
    assertEqual(n, idx)
    
    %check that error thrown when passed a list including invalid values
    bad = @() a.dval2num([2, -13]);
    assertExceptionThrown(bad, 'rpLattice:InvalidValue');
end

%% Test Dsim method
function testDsim(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice %#ok<NUSED>
    
    % call the sim helper to establish the "correct" answers
    % This helper also sets and resets the random number stream for
    % consistancy.
    seed = 4793;
    [v_test, n_test] = SimHelper(0, TLtR_tmax, TLtR_start, seed);  
        
    % Test proper handling of individual t values (t>0) including longer
    % than value_list
    for t = 0:TLtR_tmax
        [v, n] = a.dsim(t);
        assertEqual(v, v_test(t+1));
        assertEqual(n, n_test(t+1));
        assertEqual(a.t, t)
        assertEqual(a.curState(), v_test(t+1))
    end

    % Test invalid values of t (<0)
    bad = [-1, -5];
    [v, n] = a.dsim(bad);
    assertEqual(v, [NaN NaN]);
    assertEqual(n, [NaN NaN]);
    %t should keep its old value (from loop)
    assertEqual(a.t, t)
    assertEqual(a.curState(), v_test(t(end)+1))

    %--Test mixed valid & invalid values for t
    % Find "right" random based answers
    seed = 2009;
    [v_test, n_test] = SimHelper(0, TLtR_tmax, TLtR_start, seed);  

    t = 3;
    bad = [-1 2 t bad 5];
    [v, n] = a.dsim(bad);
    assertEqual(v, [NaN v_test([2 t]+1)' NaN NaN v_test(end)]);
    assertEqual(n, [NaN n_test([2 t]+1)' NaN NaN n_test(end)]);
    %t should be set to our last good value
    assertEqual(a.t, 5)
    assertEqual(a.curState(), v_test(end))


    %-- test call with vector of good values
    % Find "right" random based answers
    seed = 1974;
    [v_test, n_test] = SimHelper(0, TLtR_tmax, TLtR_start, seed);  

    t = [0 1 3];
    [v, n] = a.dsim(t);
    assertEqual(v, v_test(t+1)')
    assertEqual(n, n_test(t+1)')
    assertEqual(a.t, t(end))
    assertEqual(a.curState(), v_test(t(end)+1))

    %-- verify things work when getting values only
    % Also test a column time vector
    % Find "right" random based answers
    seed = 1492;
    v_test = SimHelper(0, TLtR_tmax, TLtR_start, seed);  

    t = [1 3 0 2 1]';
    v = a.dsim(t);
    assertEqual(v, v_test(t+1))
    assertEqual(a.t, t(end))
    assertEqual(a.curState(), v_test(t(end)+1))
    
    %--test with initial value set
    % Find "right" random based answers
    init = 15/2;
    
    seed = now;
    [v_test, n_test] = SimHelper(1, TLtR_tmax, init, seed); 

    t = [1 2 3];
    [v, n] = a.dsim(t, init);
    assertEqual(v, v_test(t)')
    assertEqual(n, n_test(t)')
    assertEqual(a.t, t(end))
    assertEqual(a.curState(), v_test(t(end)))
    
    %-- check for error with invalid value/state for specified start time
    bad = @() a.dsim([1 2 3], 20);
    assertExceptionThrown(bad, 'rpLattice:InvalidValueAtTime')
    
    
end

%% Test Sim method
function testSim(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice %#ok<NUSED>
    
    % call the sim helper to establish the "correct" answers
    % This helper also sets and resets the random number stream for
    % consistancy.
    [v_test, n_test] = SimHelper(0, TLtR_tmax, TLtR_start);  
        
    % Test proper handling of individual t values (t>0) including longer
    % than value_list
    for t = 0:TLtR_tmax
        [v, n] = a.sim(t);
        assertEqual(v, v_test(t+1));
        assertEqual(n, n_test(t+1));
        assertEqual(a.t, t)
        assertEqual(a.curState(), v_test(t+1))
    end

    % Test invalid values of t (<0)
    bad = [-10, -1, -0.2];
    [v, n] = a.sim(bad);
    assertEqual(v, [NaN NaN NaN]);
    assertEqual(n, [NaN NaN NaN]);
    %t should keep its old value (from loop)
    assertEqual(a.t, t)

    %--Test mixed valid & invalid values for t
    % Find "right" random based answers
    [v_test, n_test] = SimHelper(0, TLtR_tmax, TLtR_start);  

    t = 3;
    bad = [-1 2.1 t bad];
    [v, n] = a.sim(bad);
    assertEqual(v, [NaN v_test([2 t]+1)' NaN NaN NaN]);
    assertEqual(n, [NaN n_test([2 t]+1)' NaN NaN NaN]);
    %t should be set to our last good value
    assertEqual(a.t, t)
    assertEqual(a.curState(), v_test(t(end)+1))

    %--test call with vector of good values
    % Find "right" random based answers
    [v_test, n_test] = SimHelper(0, TLtR_tmax, TLtR_start);  

    t = [1.2 1.6 2.1 3.2 1.3]';
    [v, n] = a.sim(t);
    %in this case we want zero-order-hold & therefor use floor for
    %comparisons
    t = floor(t);
    assertEqual(v, v_test(t+1));
    assertEqual(n, n_test(t+1));
    assertEqual(a.t, t(end))
    assertEqual(a.curState(), v_test(t(end)+1))

    %-- verify things work when getting values only
    % Find "right" random based answers
    v_test = SimHelper(0, 10, TLtR_start);  

    t = [2 0.2 3.6 1.8 0.6 2.9 10];
    v = a.sim(t);
    %in this case we want zero-order-hold & therefor use floor for
    %comparisons
    t = floor(t);
    assertEqual(v, v_test(t+1)');
    assertEqual(a.t, t(end))
    assertEqual(a.curState(), v_test(end))

    %--test with initial value set
    % Find "right" random based answers
    init = 15/2;
    
    [v_test, n_test] = SimHelper(1, TLtR_tmax, init); 

    t = [1.4 2.3 3.1];
    [v, n] = a.sim(t, init);

    t = floor(t);
    assertEqual(v, v_test(t)')
    assertEqual(n, n_test(t)')
    assertEqual(a.t, t(end))
    assertEqual(a.curState(), v_test(t(end)))
    
    %-- check for error with invalid value/state for specified start time
    bad = @() a.sim([1.1 2.3 3], 20);
    assertExceptionThrown(bad, 'rpLattice:InvalidValueAtTime')
    
end

%% Test Range method
function testRange(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice

    % Test extraction of all possible states
    [v, n] = a.range('all');
    assertEqual(v, TLtR_value_list([1, end])')
    assertEqual(n, [1, length(TLtR_value_list)])
    
    % Test proper handling of t<0
    bad = @() a.range(-0.1);
    assertExceptionThrown(bad, 'RandProcess:InvalidTime');

    % Test proper handling of 0<t<Tmax
    for t = 0:TLtR_tmax
        [v, n] = a.range(t);
        assertEqual(v, [TLtR_lattice{t+1}(1) TLtR_lattice{t+1}(end)]);
        assertEqual(n, [TLtR_num_lattice{t+1}(1) TLtR_num_lattice{t+1}(end)]);
    end
    
    % Test proper handling of t>Tmax
    t_max = TLtR_tmax;
    for t = t_max+1:t_max+10
        [v, n] = a.range(t);
        assertEqual(v, [TLtR_lattice{t_max+1}(1) TLtR_lattice{t_max+1}(end)]);
        assertEqual(n, [TLtR_num_lattice{t_max+1}(1) TLtR_num_lattice{t_max+1}(end)]);
    end
    
    % Test use of current (unspecified) time
    t=2;
    
    a.t = t;
    [v, n] = a.range();
    assertEqual(v, [TLtR_lattice{t+1}(1) TLtR_lattice{t+1}(end)]);
    assertEqual(n, [TLtR_num_lattice{t+1}(1) TLtR_num_lattice{t+1}(end)]);
        
    % Test proper handling of non-integer times
    for t = [1.2 0.1 3.1 pi]
        [v, n] = a.range(t);
        t = floor(t); %#ok<FXSET>
        assertEqual(v, [TLtR_lattice{t+1}(1) TLtR_lattice{t+1}(end)]);
        assertEqual(n, [TLtR_num_lattice{t+1}(1) TLtR_num_lattice{t+1}(end)]);
    end
    
    % Test reduced number of outputs
    t = 3;
    % value only
    v = a.range(t);
    assertEqual(v, [TLtR_lattice{t+1}(1) TLtR_lattice{t+1}(end)]);
    
end

%% Test Step method
function testStep(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice %#ok<NUSED>
        
    %-- Test proper action when steping through the range of t
    % Note: the internal t starts at 0

    t_max = TLtR_tmax + 10;
    % Find "right" random based answers, including t>Tmax
    [v_test, n_test] = SimHelper(0, t_max, TLtR_start);  

    for t = 1:t_max
        [v, n, t_out] = a.step();
        assertEqual(v, v_test(t+1));
        assertEqual(n, n_test(t+1));
        assertEqual(a.t, t)
        assertEqual(t, t_out)
        assertEqual(a.curState(), v_test(t+1))
    end

    
    %-- Now go backwards
    % Note: here we keep two parallel random streams, the default one is
    % used by step and the second one we use when drawing random numbers
    % Since they are both initialized to the same value, we should get
    % identical results
    
    seed = second(now)*1000;
    rs=RandStream('mt19937ar', 'Seed', seed);
    RandStream.setGlobalStream(rs);
    rs2=RandStream('mt19937ar', 'Seed', seed);

    for t = (t_max-1):-1:0

        % Do our own backward step using listprev and the second random
        % number stream
        [v_list, s_list, prob] = a.dlistprev(); %#ok<ASGLU>
        trans = cumsum(prob);
        idx = find(rand(rs2) <= trans, 1, 'first');
        v_test = v_list(idx);

        % And compare this to the backward step
        [v, n, t_out] = a.step(-1);
        
        assertEqual(v, v_test)
        assertEqual(n, a.dval2num(v_test))
        assertEqual(a.t, t)
        assertEqual(t_out, t)
    end
    
    
    % Test backing into invalid values of t (<0)
    a.t = 0;
    for st  = [-1, -3]
        [v, n] = a.step(st);
        assertEqual(v, []);
        assertEqual(n, []);
        %t should stop at 0
        assertEqual(a.t, 0)
    end

    %-- Test bigger steps
    t = 0;

    for st = 1:3
        % Find "right" random based answers
        [v_test, n_test] = SimHelper(t, TLtR_tmax, TLtR_lattice{t+1}(t+1));  

        a.t=t;
        idx = st+1;
        [v, n, t_out] = a.step(st);
        assertEqual(v, v_test(idx));
        assertEqual(n, n_test(idx));
        assertEqual(a.t, t+st)
        assertEqual(t+st, t_out)
    end

    % Test zero stepsize (should have no change)
    t = a.t;
    st = 0;
    [v, n, t_out] = a.step(st);
    assertEqual(v, v_test(idx));
    assertEqual(n, n_test(idx));
    assertEqual(a.t, t+st)
    assertEqual(t+st, t_out)
    
    %test non-integer steps
    t = 0;
    for st = [1.2 0.1 2 pi-1]
        % Find "right" random based answers
        [v_test, n_test] = SimHelper(t, TLtR_tmax, TLtR_lattice{t+1}(t+1));  

        a.t=t;
        idx = floor(st)+1;
        [v, n, t_out] = a.step(st);
        assertEqual(v, v_test(idx));
        assertEqual(n, n_test(idx));
        assertEqual(a.t, t+st)
        assertEqual(t+st, t_out)
    end

    % Test reduced number of outputs & starting with a non-zero start time
    t = 1;
    % Find "right" random based answers
    [v_test, n_test] = SimHelper(t, TLtR_tmax, TLtR_lattice{t+1}(t+1));  

    a.t = t;
    % value only
    v = a.step();
    assertEqual(v, v_test(2));
    % value & state only
    [v, n] = a.step();
    assertEqual(v, v_test(3));
    assertEqual(n, n_test(3));

    
end

%% Test CurState method
function testCurState(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice %#ok<NUSED>

    % Test proper handling of 0<t<Tmax
    for t = 0:(length(TLtR_lattice)-1)
        a.t = t;
        [v, n, t_out] = a.curState();
        idx = ceil(length(TLtR_lattice{t+1})/2);
        assertEqual(v, TLtR_lattice{t+1}(idx));
        assertEqual(n, TLtR_num_lattice{t+1}(idx));
        assertEqual(t, t_out);
    end
    
    % Test reduced number of outputs
    t = 2;
    a.t = t;
    % value only
    v = a.curState();
    idx = ceil(length(TLtR_lattice{t+1})/2);
    assertEqual(v, TLtR_lattice{t+1}(idx));
    % value & state only
    [v, n] = a.curState();
    assertEqual(v, TLtR_lattice{t+1}(idx));
    assertEqual(n, TLtR_num_lattice{t+1}(idx));
end

%% Test Reset method
function testReset(a)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice %#ok<NUSED>

    a.reset();
    
    assertEqual(a.t, 0);
    assertEqual(a.curState(), TLtR_start);
end





%% ---------------------------- Helper Functions -------------------------------
function CheckPrevStates(a, t, num_out)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice

    if nargin < 2
        t = a.t;
        specify_t = false;
    else
        specify_t = true;
    end
    
    if nargin < 3
        num_out = 3;
    end
        
    
    % For all possible states in this time period
    for s_idx = 1:length(TLtR_num_lattice{t+1})
        s=TLtR_num_lattice{t+1}(s_idx);
        this_v = TLtR_value_list(s);
        if specify_t
            if num_out == 3
                [prev_v, n, p] = a.dlistprev(s,t);
            elseif num_out == 2
                [prev_v, n] = a.dlistprev(s,t);
            else
                prev_v = a.dlistprev(s,t);
            end
        else
            if num_out == 3
                [prev_v, n, p] = a.dlistprev(s);
            elseif num_out == 2
                [prev_v, n] = a.dlistprev(s);
            else
                prev_v = a.dlistprev(s);
            end
        end

        % make sure that all possible prior states are listed with the
        % correct probabilities by looping over all transition
        % coeficients
        vals_found = 0;
        for c_idx = 1:length(TLtR_coef)
            possible_v = RoundTo(this_v/TLtR_coef(c_idx), 0.0001);
            [found_v, found_idx] = intersect(possible_v, TLtR_lattice{t}); %#ok<NASGU>
            if isempty(found_v)
                continue
            elseif length(found_v) > 1
                %If we get here, our lattice is wrong
                error('testrpLattice:toomanyprev', 'Multiple previous states found for one coefficient')
            else
                vals_found = vals_found+1;
                assertTrue(ismember(found_v, prev_v))
%TODO: Check for proper probabilities
            end
            assertTrue(vals_found>0)
        end
        if num_out >= 2
            assertEqual(TLtR_value_list(n), prev_v);
        end
        if num_out == 3
            assertEqual(sum(p), 1)
        end
    end
end

function [v_test, n_test] = SimHelper(t_min, t_max, cur_val, seed)
    global TLtR_start TLtR_coef TLtR_prob TLtR_tmax TLtR_lattice %#ok<NUSED>
    global TLtR_value_list TLtR_num_lattice 

    if nargin < 4 || isempty(seed)
        seed = second(now)*1000;
    end
    
    % reset the random stream using a seed-based rand sequence
    % Note: we assume the rands are called in simulation order and that all
    % rand values are used. This may have to change with alternate
    % implementations
    rs=RandStream('mt19937ar', 'Seed', seed);
    RandStream.setGlobalStream(rs);
    
    v_test = NaN(t_max-t_min+1, 1);
    n_test = NaN(size(v_test));
    v_test(1) = cur_val;
    n_test(1) = find(v_test(1) == TLtR_value_list);
    
    for t_step = 2: (t_max-t_min+1)
        t_lookup = min(t_step, TLtR_tmax+1);
        if t_step > TLtR_tmax+1
            v_test(t_step) = v_test(t_step-1);
        else
            next = find(rand <= cumsum(TLtR_prob), 1, 'first');
            v_test(t_step) = RoundTo(v_test(t_step-1)*TLtR_coef(next), 0.0001);
        end
        n_test(t_step) = find(v_test(t_step) == TLtR_lattice{t_lookup+t_min});
        n_test(t_step) = TLtR_num_lattice{t_lookup+t_min}(n_test(t_step));
    end

    % reset the random stream using a seed-based rand sequence
    rs=RandStream('mt19937ar', 'Seed', seed);
    RandStream.setGlobalStream(rs);

end
