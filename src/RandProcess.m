classdef RandProcess < AbstractSet
%RANDPROCESS Random Process abstract class for dynamic programming
%
% Defines an abstract class (defines the strucuture for related
% subclasses) for a Random Process for use with dynamic programming and
% other simulations
%
% IMPORTANT: all processes must have a single, starting state for t=0
%
% originally by Bryan Palmintier 2010

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   8  2016-10-07 01:40  BryanP      Convert as_array to wrapper for dlistnext 
%   7  2016-10-06 11:40  DheepakK    as_array that throws an error
%   6  2016-07-07 09:40  BryanP      Extracted out common pieces to form AbstractSet class
%   5  2010-12-23 19:30  BryanP      Made t an abstract property to allow subclass error handling
%   4  2010-12-15 23:30  BryanP      Added t_max for dlist(... 'all')
%   3  2010-12-14 12:44  BryanP      Added value returns to dlist*
%   2  2010-12-13 21:20  BryanP      distinguish dsim & sim add internal t
%   1  2010-12-13 12:20  BryanP      Initial Version

    properties (Abstract)
        t  %current timestep
    end

    methods (Abstract)
        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.
        %
        % IMPORTANT: all processes must have a single, starting state for t=0
        % DLIST List possible discrete states
        %
        % List possible discrete states by number for given time
        % if t is not listed, the states for the current simulation time
        % are returned.
        %
        % To get a list of all possible states pass with t='all', in which
        % case, t_max is provided to bound the problem size.
        [value_list, state_n_list] = dlist (obj, t)

        % DLISTPREV List previous discrete states & probabilities
        %
        % List possible previous states (by number) along with conditional
        % probability P(s_t | s_{t-1})
        %
        % If t is not defined, the current simulation time is assumed
        [value_list, state_n_list, prob] = dlistprev (obj, state_n, t )

        % DLISTNEXT List next discrete states & probabilities
        %
        % List possible next states (by number) along with conditional
        % probability P(s_{t+1} | s_t)
        %
        % If t is not provided, the current simulation time is assumed
        %
        % If the state is not valid at t, an error with ID
        % 'RandProcess:InvalidState' should be thrown
        %
        % If t is out of the valid range, an error with ID
        % 'RandProcess:InvalidTime'
        [value_list, state_n_list, prob] = dlistnext (obj, state_n, t )

        % DNUM2VAL Convert discrete state number to a value
        values = dnum2val (obj, state_n_list )

        % DVAL2NUM Convert the value(s) to associate discrete state numbers
        state_n_list = dval2num (obj, values )

        % DSIM Simulate discrete process.
        %
        % A column vector for t is assumed to be a series of times for
        % which to return results. Intermediate times are also computed, if
        % needed, but not returned. The initial value is not returned in
        % the value series. Only one simulation is run, such that out of
        % order times will be sorted before simulation and duplicate times
        % will return the same result
        %
        % Invalid times (t<=0) return NaN
        %
        % Note: after calling sim, the process internal time will be set to
        % the final value of t_list
        [value_series, state_n_series] = dsim(obj,t_list, initial_value)

        %% ===== General (discrete or continuous) Methods

        % SIM Simulate process for desired (continuous) times
        %
        % A column vector for t is assumed to be a series of times for
        % which to return results. Intermediate times are also computed, if
        % needed. The initial value is not returned in the value series.
        %
        % Function must handle arbitrary positive values for t_list
        % Invalid times return NaN
        %
        % Note: after calling sim, the process internal time will be set to
        % the final value of t_list
        [value_series, state_n_series] = sim(obj, t_list, initial_value)

        % RANGE Find value range for given time
        %
        % Returns vector with [min max] value range for specified time
        % if t is not provided, the range for the current simulation time
        % is returned.
        %
        % To get the possible range across all times use t='all'
        [value_range, state_n_range] = range(obj, t)

        %% ===== Additional simulation support
        %STEP simulate forward
        %
        % by default steps forward by delta_t = 1
        [value, state_n, t] = step(obj, delta_t)

        %CURSTATE Return the current state of the simulation
        [value, state_n, t] = curState(obj)

        %RESET reset simulation to t=0 and provided intitial_value
        reset(obj, initital_value)

    end

    methods

        function [val, prob] = as_array(obj)
            %For a RandProcess, as_array is a simple wrapper around
            %dlistnext for the current state and time.c

            [val, ~, prob] = obj.dlistnext(obj.cur_state, obj.t);

        end

    end

end
