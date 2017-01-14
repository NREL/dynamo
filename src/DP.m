function [policy, values, params, states] = DP(n_periods, disc_rate, fStateIterator, ...
            fDecisionIterator, fFinalValue, fContribution, fTransProbs, ...
            params, verbose)
% DP Solve Dynamic Program (ESD.862)
%
% Solves generic dynamic programming problems by backward induction,
% assuming the goal is to maximize the discounted expected value
%
% Note: backward induction assumes a markov relation between states
%
% Usage: [policy, values, params, states] = DP(n_periods, disc_rate, @fStateIterator, ...
%                @fDecisionIterator, @fFinalValue, @fContribution, @fTransProbs,...
%                params,verbose)
%
% PARAMETERS... (all required)
%   n_periods  the number of time periods
%   disc_rate  the discount rate
%
%   Functions & Parameters
%    the core of DP call is in specifyin problem specific functions for
%    each step of the algorithm. Specifically the user must provide
%    functions of the following forms:
%
%     [s, params] = fStateIterator(i, t, params)
%           returns a state variable (of any type), for the positive
%           decimal state index "i". Specifically:
%              - i = 0: returns the max number of possible states
%              - 0 < 1 < number of states: returns a state variable of any
%                valid Matlab type. If the return is NaN, the state is
%                assumed to be invalid (at the current t) and processing
%                skips to the next state
%     [d, params] = fDecisionIterator(s, i, t, future_values, params)
%           returns a decision, represented by a variable of any type, that
%           represents one of the possible decisions to be made at state s.
%                i = '0 returns the number of possible decisions for the
%                current state. Otherwize the range of i is 1:num decisions
%                for this state.
%           The future_values parameter allows drastically reduced decision
%           spaces based on unreachable, or otherwise undesirable next
%           states.
%     [v, params] = fFinalValue(s, t, params)
%           returns the final value for state s
%     [c, params] = fContribution(s, d, t, params)
%           the contribution this period of making decision d in given
%           state, s at time/solution step t. In the inventory problem,
%           this combines ordering & holding costs. Note costs should be
%           negative
%     [p_vec, c_vec, params] = fTransProbs(s, d, t, params)
%           Computes the transition probability and contributions for the
%           decision d in timestep (or solution step) t and returns 2 column
%           vectors:
%        p_vec = a vector of probabilies such that p_vec(i) corresponds
%                to the probability of transitioning to state s+1, given by
%                state_iterator(i) from state s for the given decision.
%                  length(p_vec) = state_iterator('end')
%        c_vec = the cost (use negative for gain) associated with the
%                corresponding transition. In the inventory problem this is
%                the negative income from sales.
%     OR
%     [trans_struct, ~, params] = fTransProbs(s, d, t, params)
%           Where trans_struct contains the fields:
%        trans_struct.prob with all of the non-zero probabilities in p_vec
%        trans_struct.cost with corresponding non-zero costs in c_vec
%        trans_struct.idx  a vector of corresponding linear state indices
%
%
%    params    additional parameter to be used by the called functions.
%              Note: use a struct to pass multiple parameters
%  OPTIONAL PARAMS...
%    verbose   a scalar indicating the level of verbose output, a value of
%              zero implies no output. The number indicates which itertion
%              multiples to update the display (default = 0)
%
%  RETURNS...
%   policy    a cell array of all optimal decisions. For simple numeric
%             decisions, convert back to a standard array with cell2mat
%   values    a numeric array of the optimal value of each state in each
%             time period.
%   states    a cell array showing all of the corresponding states for the
%             policy and values.
%
% Originally by Bryan Palmintier, 2010


% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-10 15:47  BryanP      initial adaptation from inventory_dp version 3
%   2  2010-05-11 21:17  BryanP      completed adaptation & initial tests
%   3  2010-05-14 20:00  BryanP      converted to external param intialization
%   4  2010-05-14 20:30  BryanP      now i==0 for number of states/decisions
%   5  2010-05-22 09:50  BryanP      added support for time-aware states
%   6  2010-05-23 00:40  BryanP      time-aware decisions & pass terminal
%                                     state number to finalValue function
%   7  2010-06-05 20:15  BryanP      Store NaN for invalid state policy & value
%   8  2010-06-06 01:00  BryanP      Optional verbose status display
%   9  2010-06-07 00:30  BryanP      Allow param update with each sub-function call
%  10  2010-06-07 08:00  BryanP      Pass future state vector to fDecision
%  11  2010-06-07 08:00  BryanP      Pass future state vector to fTransProb
%  12  2010-06-08 01:26  BryanP      Corrected/faster expected value for NaN
%  13  2011-03-10 09:05  BryanP      Abort early if no valid states in time period
%  14  2011-04-17 23:05  BryanP      Support for trans prob with explict indices rather than full vectors
%  15  2011-06-08 16:15  BryanP      Added isscalar() to NaN state checks: avoid err with numeric vector states
%  16  2011-07-11 10:00  BryanP      abstracted out DisplayProgress for '.' display
% ----- Handle optional inputs -----
if nargin < 9 || isempty(verbose)
    verbose = 0;
end

% ----- Setup Problem -----
% -- Parameter intialization
%Determine the state vector size
[n_states, params] = fStateIterator(0, [], params);

if verbose
    fprintf('\nDP: setting up. %d possible states\n',n_states)
end

%Note: we are keeping track of the value matrix right here in this
%function, so we need to initialize it...

%value function = f(s,t). Initialize to -Inf so any real value is greater
values = -Inf * ones(n_states, n_periods+1);

%initialize policy array
%Note: as with the value array, we (this function) is in charge of this
%matrix. But unlike the value matrix, we want to allow the user to use any
%type of data, not just a single number, to describe the decision,
%therefore we will use a cell array.
policy = NaN * ones(n_states, n_periods); %action (orders) = f(s,t), no action in final state
policy = num2cell(policy);

%-- Handle optional outputs

%If the user wants a list of states, create a flag for future reference
list_states = nargout > 2;
% and initialize the state list as another cell arrray so we can handle
% non-numeric states (like the decision/policy array)
if list_states
    states = NaN * ones(n_states, n_periods);
    states = num2cell(states);
end

%% -- compute terminal values.
if verbose
    fprintf('  Terminal Values')
end

for s_idx = 1: n_states;
    DisplayProgress(verbose, s_idx, [])

    [s, params] = fStateIterator(s_idx,n_periods+1,params);
    if isnumeric(s) && isscalar(s) && isnan(s)
        values(s_idx, n_periods+1) = NaN;
    else
        [values(s_idx, n_periods+1), params] = ...
            fFinalValue(s, n_periods+1, params);
    end
end

valid_states = nnz(not(isnan(values(:,n_periods+1))));

if verbose
    fprintf('Done: %d valid terminal states\n', valid_states)
end

if valid_states == 0
    error('DP:NoValidStates', 'No valid states for terminal period: Aborting')
end

%% --- Work backwards ---
% Note: terminal values already computed above

% Loop over time periods in reverse (backward induction)
for t = n_periods:-1:1
    if verbose
        fprintf('  Period #%d', t);
    end
    % Loop over states
	for s_idx = 1:n_states;
        DisplayProgress(verbose, s_idx, [])

		[s, params] =fStateIterator(s_idx, t, params);
        if list_states
            states{s_idx, t} = s;
        end
        %Skip processing if state is invalid (used to skip unreachable
        %states at certain times
        if isnumeric(s) && isscalar(s) && isnan(s)
            values(s_idx, t) = NaN;
            policy{s_idx, t} = NaN;
            continue
        end

        %Now Loop over each possible decision for this state
        % First find the number of decisions by passing decision index = 0
        [n_decisions, params] = fDecisionIterator(s, 0, t, values(:,t+1), params);

        % Now loop over all decisions
        for d_idx = 1:n_decisions;
            [d, params] = fDecisionIterator(s, d_idx, t, values(:,t+1), params);

            % compute immediate decision contribution
            [contrib_now, params] = fContribution(s, d, t, params);

            % find transition probabilities and corresponding near-term
            % costs
            [p_vec, c_vec, params] = fTransProbs(s, d, t, values(:,t+1), params);

            % compute expected value of this decisions as
            %  imediate cost + expected near-term costs + expected value of
            %  next state

            % NEW style transitions (pass only non-zero elements and corresponding indicies):
            if isstruct(p_vec)
                exp_value = contrib_now + p_vec.prob' ...
                              * (p_vec.cost + (1-disc_rate) .* values(p_vec.idx,t+1));

            % OLD style transitions (pass one p & c element for each state):
            else
                %extract indicies of valid next states
                possible = find(p_vec > 0);

                %compute expected value from valid states
                exp_value = contrib_now + p_vec(possible)' ...
                              * (c_vec(possible) + (1-disc_rate) .* values(possible,t+1));
            end

            if  exp_value > values(s_idx, t)
				values(s_idx, t) = exp_value;
				policy{s_idx, t} = d;
            end
        end %decision

	end 	%End (loop over all states s for time t)

    valid_states = nnz(not(isnan(values(:,t))));

    if verbose
        fprintf('Done: %d valid states\n', valid_states)
    end

    if valid_states == 0
        error('DP:NoValidStates', 'No valid states for period %d: Aborting', t)
    end

end % End (loop over all time periods)
