classdef rpTransMatrix < RandProcess
%rpTransMatrix arbitrary n-to-m-to-... transition defined by set of matrices
%
% Generic random process with fully specified n by m transition matricies
% per time.
%
% Constructor: obj = rpTransMatrix(states, trans, t_max)
%
% Required inputs/properties:
%  states: cell row vector with a set of states for each time period.
%    Each state set is a column-wise list where each row represents a
%    state.
%  trans:  cell row vector with a set of transition matrices for each time
%    period. Each matrix should be of size n_states(t) x n_states(t+1).
% Optional inputs/properties:
%  t_max: maximum allowable time. Defaults to length(states). If t_max >
%    length(states), the final transition matrix (trans{end}) must be
%    square
%
% Notes:
%  - Integer timesteps are used, scale accordingly for fractional
%    timesteps
%  - start is the value at t=1
%  - For non-integer times, the output is assumed to remain at the value
%    corresponding to the previous integer time step (Zero-order hold)
%  - Times t>t_max (only valid for square final transition matricies)
%    become Markov processes with a fixed set of states 
%
% Examples (using ISGT paper values)
%
% >> rng(0);
% >> states = {0.76, [0.76+0.33; 0.76+0.1], [0.76+0.33+0.42; .76+0.1+0.4; 0.76+0.1+0.25]};
% >> prob = {[1 1]/2, [0.80 0.20  0   ; 0   0.2 0.8]};
% >> rtm_obj = rpTransMatrix(states, prob)
% 
% rtm_obj = 
% 
%   rpTransMatrix with properties:
% 
%              Tol: 1.0000e-06
%                t: 1
%        cur_state: 0.7600
%             name: ''
%     pt_dim_names: {}
%       SampleType: 'rand'
%            N_dim: 1
%     DiscreteMask: []
%
% >>[s_list, p_list] = rtm_obj.as_array(1)
% s_list =
%     0.7600
% p_list =
%      1
%
% >> [s_list, p_list] = rtm_obj.as_array(3)
% s_list =
%     1.5100
%     1.2600
%     1.1100
% p_list =
%     0.4000
%     0.2000
%     0.4000
% 
% >> rtm_obj.sample(3)
% ans =
%     0.7600
%     0.7600
%     0.7600
%
% >> rtm_obj.sample(5, 3)
% ans =
%     1.1100
%     1.1100
%     1.5100
%     1.5100
%     1.2600
% 
%
% see also rpLattice, RandProc, rpDiscreteSample, rpMarkov, rpBasic
%
% originally by Bryan Palmintier 2017

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   2  2017-07-16 00:13  BryanP      Use standardized conditionatlSample() in RandProcess 
%   1  2017-07-15 20:02  BryanP      Original version adapted from rpLattice v15

    % Internal properties
    properties (Access='protected')
        Transitions = {};   % set of transition matrices
    end

    methods
        %% ===== Constructor & related
        function obj = rpTransMatrix(states, trans, t_max)
            % Note: see rpLattice class documentation (above) for more info
            %
            % Support zero parameter calls to constructor for special
            % MATLAB situations (see help files)
            if nargin == 0
                return
            end
            if nargin < 3 || isempty(t_max)
                t_max = length(states);
            end
            
            %-- Check for consistant inputs
            % Should have one more state_list than transition probability
            if length(states) > length(trans) + 1
                warning('RandProc:LengthMismatch', 'Too many states (%d), truncating to one more than number of transition matrices (%d)', length(states), length(trans) + 1);
                states( (length(trans)+1):end ) = [];
            elseif length(trans) > length(states) - 1
                warning('RandProc:LengthMismatch', 'Too many transition matrices (%d), truncating to one less than number of states (%d)', length(trans), length(states) - 1);
                trans( (length(states)-1):end ) = [];
            end
            
            % Allow only one initial state
            if size(states{1}, 1) ~= 1
                warning('RandProc:NonDeterministicStart', 'More than one intial state defined, keeping only the first ([ %s])', sprintf('%g ',states{1}));
                states{1}( 2:end ) = [];
            end
            
            % Verify trans dimensions match states & probability sum to one
            % (by row)
            %  Each matrix should be of size n_states(t) x n_states(t+1)
            %  Note: we will transpose before multiplying
            for t = 1:length(trans)
                % Check dimensions
                assert( isequal( size(trans{t}), [length(states{t}), length(states{t+1})] ), ...
                    'RandProc:InvalidSize', ...
                    'transition matrix size ([ %s]) must equal size n_states(t) x n_states(t+1) or ([ %d %d ]) at t=%d', ...
                    sprintf('%g ', size(trans{t})), length(states{t}), length(states{t+1}) )
                
                % Check probabilites sum to one (by row)
                assert( not(any( abs(sum(trans{t}, 2) - 1) > obj.Tol)), ...
                    'RandProc:ProbSum1', ...
                    'Transition matrix probabilities must sum to 1.0 by row for t=%d', t )
                
            end
            
            % Support t > length(states) only if final transition is square
            if t_max > length(states)
                [nrow, ncol] = size(trans{end});
                if nrow ~= ncol
                    t_max = length(states);
                    warning('rpTransMatrix:TmaxMismatch', 'Limiting tmax to last t with defined states (%d) because final transition not square', t_max);
                end
            end
            
            %--- Store parameters
            obj.Values = states;
            obj.Transitions = trans;
            obj.Tmax = t_max;
            
            % Compute and store unconditional probabilities
            obj.UncondProbs{1} = 1;
            obj.UncondCdfs{1} = 1;
            for t = 2:obj.N_uniqueT
                obj.UncondProbs{t} = obj.UncondProbs{t-1}' * obj.Transitions{t-1};
                obj.UncondProbs{t} = obj.UncondProbs{t}';
                
                obj.UncondCdfs{t} = cumsum(obj.UncondProbs{t});
            end
            
            obj.reset();
        end
        
        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.

        function [state_list, prob] = dlistnext (obj, state_in, t )
        % DLISTNEXT List next discrete states & probabilities
        %
        % List possible next states (by number) along with conditional
        % probability P(s_{t+1} | s_t)
        %
        % If t and/or state are not defined, the current simulation time
        % and state are assumed
            if nargin < 2 || isempty(state_in)
                state_in = obj.cur_state;
            end
            
            if nargin < 3 || isempty(t)
                t = obj.t;
            end
            
            [state_idx, t_lookup] = obj.checkState(t, state_in);
            
            % Find next value list
            next_prob = obj.Trans{t_lookup}(state_idx, :)';
            valid_state_map = next_prob > 0;
            state_list = obj.Values{t_lookup}(valid_state_map);

            if nargout > 1
                prob = next_prob(valid_state_map);
            end
        end
    end
    
    methods (Access = protected)
    end

end