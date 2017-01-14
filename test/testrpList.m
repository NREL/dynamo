function test_suite = testrpList %#ok<MCUOA>
%TESTrpList Test functions for rpList class
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
%   1  2010-12-13 23:30  BryanP      Original Code
%   2  2010-12-14        BryanP      Expanded to cover more functions
%   3  2010-12-15 23:30  BryanP      Adapted for list start at t=0
%   4  2011-01-04 22:00  BryanP      Minor tweaks to match rpList v9

    initTestSuite;
end

%% Shared Setup, used by all tests
function a = setup %#ok<*DEFNU>
    global TLR_value_list TLR_state_num_list TLR_possible_value_list
    global TLR_value_series TLR_state_num_series

    TLR_value_list =    [0 2 2 10 3 2 10 15];
    TLR_state_num_list = [1 2 2  4 3 2  4 5];
    TLR_possible_value_list = [0 2 3 10 15];
    TLR_value_series =     [TLR_value_list     15 15 15 15];
    TLR_state_num_series = [TLR_state_num_list  5  5  5  5];
    
    a = rpList(TLR_value_list);
end

%% Shared Teardown, used by all tests
function teardown(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series %#ok<NUSED>

    a.delete;
    clearvars a TLR_*
end

%% Test constructor
function testConstructor(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series %#ok<NUSED>

    assertEqual(a.ValueSeries, TLR_value_list)
    assertEqual(a.t, 0)
end

%% Test Dlist method
function testDlist(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series

    % Test extraction of all possible states
    [v, n] = a.dlist('all');
    assertEqual(v, TLR_possible_value_list)
    assertEqual(n, 1:length(TLR_possible_value_list))
    
    % Test proper handling of t<0
    bad = @() a.dlistprev(-1);
    assertExceptionThrown(bad, 'RandProcess:InvalidTime')
            
    % Test proper handling of t>=0, including t> length of value_list
    for t = 0:(length(TLR_value_series)-1)
        [v, n] = a.dlist(t);
        assertEqual(v, TLR_value_series(t+1));
        assertEqual(n, TLR_state_num_series(t+1));
    end
    
    % Test use of current (unspecified) time
    t=2;
    
    a.t = t;
    [v, n] = a.dlist();
    assertEqual(v, TLR_value_series(t+1));
    assertEqual(n, TLR_state_num_series(t+1));

    % Test reduced number of outputs
    t = 5;
    a.t = t;
    % value only
    v = a.dlist();
    assertEqual(v, TLR_value_series(t+1));
end

%% Test Dlistprev method
function testDlistPrev(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series

    % Test proper handling of t<0
    bad = @() a.dlistprev([],-1);
    assertExceptionThrown(bad, 'RandProcess:InvalidTime')
    
    % Test proper handling of t=0
    bad = @() a.dlistprev([],0);
    assertExceptionThrown(bad, 'RandProcess:InvalidTime')

    
    % Test proper handling of t>=1, including t> length of value_list
    for t = 1:length(TLR_value_series)
        [v, n, p] = a.dlistprev([],t);
        assertEqual(v, TLR_value_series(t));
        assertEqual(n, TLR_state_num_series(t));
        assertEqual(p, 1)
    end
    
    % Test use of current (unspecified) time & state
    t=2;
    
    a.t = t;
    [v, n, p] = a.dlistprev();
    assertEqual(v, TLR_value_series(t));
    assertEqual(n, TLR_state_num_series(t));
    assertEqual(p, 1)
    
    % Test reduced number of outputs
    t = 5;
    a.t = t;
    % value only
    v = a.dlistprev();
    assertEqual(v, TLR_value_series(t));
    % value & state only
    [v, n] = a.dlistprev();
    assertEqual(v, TLR_value_series(t));
    assertEqual(n, TLR_state_num_series(t));
end

%% Test Dlistnext method
function testDlistNext(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series

    % Test proper handling of t<-1 
    bad = @() a.dlistprev([],-2);
    assertExceptionThrown(bad, 'RandProcess:InvalidTime')
    
    % Test proper handling of t>0, including t> length of value_list
    for t = -1:length(TLR_value_series)-2
        [v, n, p] = a.dlistnext([],t);
        assertEqual(v, TLR_value_series(t+2));
        assertEqual(n, TLR_state_num_series(t+2));
        assertEqual(p, 1)
    end
    
    % Test use of current (unspecified) time & state
    t=2;
    
    a.t = t;
    [v, n, p] = a.dlistnext();
    assertEqual(v, TLR_value_series(t+2));
    assertEqual(n, TLR_state_num_series(t+2));
    assertEqual(p, 1)

    % Test reduced number of outputs
    t = 5;
    a.t = t;
    % value only
    v = a.dlistnext();
    assertEqual(v, TLR_value_series(t+2));
    % value & state only
    [v, n] = a.dlistnext();
    assertEqual(v, TLR_value_series(t+2));
    assertEqual(n, TLR_state_num_series(t+2));
    
end

%% Test Dnum2val method
function testDnum2val(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series %#ok<NUSED>
    
    %test call with single, valid state number
    for idx = 1:length(TLR_possible_value_list)
        v = a.dnum2val(idx);
        assertEqual(v, TLR_possible_value_list(idx))
    end
    
    %test call with vector of state numbers
    idx = [1 3 1 5];
    v = a.dnum2val(idx);
    assertEqual(v, TLR_possible_value_list(idx))
    
    %check that error thrown when passed a list including an invalid state num
    bad = @() a.dnum2val([2 10 15]);
    assertExceptionThrown(bad, 'rpList:InvalidStateNum');
end

%% Test Dval2num method
function testDval2num(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series %#ok<NUSED>
    
    %test call with single, valid values
    for idx = 1:length(TLR_possible_value_list)
        v = TLR_possible_value_list(idx);
        n = a.dval2num(v);
        assertEqual(n, idx)
    end
    
    %test call with vector of values
    idx = [1 4 1 5];
    v = TLR_possible_value_list(idx);
    n = a.dval2num(v);
    assertEqual(n, idx)
    
    %check that error thrown when passed a list including invalid values
    bad = @() a.dval2num([2, -13]);
    assertExceptionThrown(bad, 'rpList:InvalidValue');
end

%% Test Dsim method
function testDsim(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series
        
    % Test proper handling of individual t values (t>0) including longer
    % than value_list
    for t = 0:(length(TLR_value_series)-1)
        [v, n] = a.dsim(t);
        assertEqual(v, TLR_value_series(t+1));
        assertEqual(n, TLR_state_num_series(t+1));
        assertEqual(a.t, t)
    end

    % Test invalid values of t (<0)
    bad = [-1, -5];
    [v, n] = a.dsim(bad);
    assertEqual(v, [NaN NaN]);
    assertEqual(n, [NaN NaN]);
    %t should keep its old value (from loop)
    assertEqual(a.t, t)

    % Test mixed valid & invalid values for t
    t = 3;
    bad = [-1 2 t bad];
    [v, n] = a.dsim(bad);
    assertEqual(v, [NaN TLR_value_series([2 t]+1) NaN NaN]);
    assertEqual(n, [NaN TLR_state_num_series([2 t]+1)    NaN NaN]);
    %t should keep its old value (from loop)
    assertEqual(a.t, t)


    %test call with vector of good values
    t = [1 4 1 5];
    [v, n] = a.dsim(t);
    assertEqual(v, TLR_value_series(t+1));
    assertEqual(n, TLR_state_num_series(t+1));
    assertEqual(a.t, t(end))

    %verify things work when getting values only
    t = [2 4 3 8 4];
    v = a.dsim(t);
    assertEqual(v, TLR_value_series(t+1));
    assertEqual(a.t, t(end))

end

%% Test Sim method
function testSim(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series
        
    % Test proper handling of individual integer t values (t>1) including
    % longer than value_list
    for t = 0:(length(TLR_value_series)-1)
        [v, n] = a.sim(t);
        assertEqual(v, TLR_value_series(t+1));
        assertEqual(n, TLR_state_num_series(t+1));
        assertEqual(a.t, t)
    end

    % Test invalid values of t (<=0)
    bad = [-10, -1, -0.2];
    [v, n] = a.sim(bad);
    assertEqual(v, [NaN NaN NaN]);
    assertEqual(n, [NaN NaN NaN]);
    %t should keep its old value (from loop)
    assertEqual(a.t, t)

    % Test mixed valid & invalid values for t
    t = 3;
    bad = [-1 2 t bad];
    [v, n] = a.sim(bad);
    assertEqual(v, [NaN TLR_value_series([2 t]+1)     NaN NaN NaN]);
    assertEqual(n, [NaN TLR_state_num_series([2 t]+1) NaN NaN NaN]);
    %t should keep its old value (from loop)
    assertEqual(a.t, t)


    %test call with vector of good values
    t = [1.2 4.6 1.8 5.7 11.3];
    [v, n] = a.sim(t);
    %in this case we want zero-order-hold & therefor use floor for
    %comparisons
    t = floor(t);
    assertEqual(v, TLR_value_series(t+1));
    assertEqual(n, TLR_state_num_series(t+1));
    assertEqual(a.t, t(end))

    %verify things work when getting values only
    t = [2 4.2 3.6 8 10.6 4.9];
    v = a.sim(t);
    %in this case we want zero-order-hold & therefor use floor for
    %comparisons
    t = floor(t);
    assertEqual(v, TLR_value_series(t+1));
    assertEqual(a.t, t(end))

end

%% Test Range method
function testRange(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series

    % Test extraction of all possible states
    [v, n] = a.range('all');
    assertEqual(v, TLR_possible_value_list([1, end]))
    assertEqual(n, [1, length(TLR_possible_value_list)])
    
    % Test proper handling of t<0
    bad = @() a.range(-0.1);
    assertExceptionThrown(bad, 'RandProcess:InvalidTime');

    % Test proper handling of t>0, including t> length of value_list
    for t = 0:(length(TLR_value_series)-1)
        [v, n] = a.range(t);
        assertEqual(v, TLR_value_series([t, t]+1));
        assertEqual(n, TLR_state_num_series([t, t]+1));
    end
    
    % Test use of current (unspecified) time
    t=2;
    
    a.t = t;
    [v, n] = a.range();
    assertEqual(v, TLR_value_series([t, t]+1));
    assertEqual(n, TLR_state_num_series([t, t]+1));
        
    % Test proper handling of non-integer times
    for t = [1.2 0.1 10.1 pi]
        [v, n] = a.range(t);
        assertEqual(v, TLR_value_series(floor(t+1)).* [1 1]);
        assertEqual(n, TLR_state_num_series(floor(t+1).* [1 1]));
    end
    
    % Test reduced number of outputs
    t = 5;
    % value only
    v = a.range(t);
    assertEqual(v, TLR_value_series(t+1 .* [1 1]));
    
end

%% Test Step method
function testStep(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series
        
    % Test proper action when steping through (and beyond) the range of t
    % Note the internal t starts at 0
    for t = 1:(length(TLR_value_series)-1)
        [v, n, t_out] = a.step();
        assertEqual(v, TLR_value_series(t+1));
        assertEqual(n, TLR_state_num_series(t+1));
        assertEqual(a.t, t)
        assertEqual(t, t_out)
    end

    % Now go backwards
    for t = (length(TLR_value_series)-2):-1:0
        [v, n, t_out] = a.step(-1);
        assertEqual(v, TLR_value_series(t+1));
        assertEqual(n, TLR_state_num_series(t+1));
        assertEqual(a.t, t)
        assertEqual(t, t_out)
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

    % Test bigger steps
    t = 0;
    a.t=t;
    for st = 2:4
        t = t + st;
        [v, n, t_out] = a.step(st);
        assertEqual(v, TLR_value_series(t+1));
        assertEqual(n, TLR_state_num_series(t+1));
        assertEqual(a.t, t)
        assertEqual(t, t_out)
    end

    % Test zero stepsize
    st = 0;
    [v, n, t_out] = a.step(st);
    assertEqual(v, TLR_value_series(t+1));
    assertEqual(n, TLR_state_num_series(t+1));
    assertEqual(a.t, t)
    assertEqual(t, t_out)
    
    %test non-integer steps
    t = 0;
    a.t=t;
    for st = [1.2 0.1 2 pi]
        t = t + st;
        [v, n, t_out] = a.step(st);
        assertEqual(v, TLR_value_series(floor(t+1)));
        assertEqual(n, TLR_state_num_series(floor(t+1)));
        assertEqual(a.t, t)
        assertEqual(t, t_out)
    end

    %verify things work when getting values only
    a.t = 0;
    v = a.step();
    assertEqual(v, TLR_value_series(1+1));
    assertEqual(a.t, 1)

    % Test reduced number of outputs
    t = 5;
    a.t = t;
    % value only
    v = a.step();
    assertEqual(v, TLR_value_series(t+2));
    % value & state only
    [v, n] = a.step();
    assertEqual(v, TLR_value_series(t+3));
    assertEqual(n, TLR_state_num_series(t+3));

    
end

%% Test CurState method
function testCurState(a)
    global TLR_value_list TLR_state_num_list TLR_possible_value_list %#ok<NUSED>
    global TLR_value_series TLR_state_num_series

    % Test proper handling of t>=0, including t> length of value_list
    for t = 0:(length(TLR_value_series)-1)
        a.t = t;
        [v, n, t_out] = a.curState();
        assertEqual(v, TLR_value_series(t+1));
        assertEqual(n, TLR_state_num_series(t+1));
        assertEqual(t, t_out);
    end
    
    % Test reduced number of outputs
    t = 5;
    a.t = t;
    % value only
    v = a.curState();
    assertEqual(v, TLR_value_series(t+1));
    % value & state only
    [v, n] = a.curState();
    assertEqual(v, TLR_value_series(t+1));
    assertEqual(n, TLR_state_num_series(t+1));
end

%% Test Reset method
function testReset(a)
    a.reset();
    
    assertEqual(a.t, 0);
end
