classdef faDiscrete < FuncApprox
%Updated for NEW FuncApprox starting from April 17th 2016
%
%
%FADISCRETE Discrete function approximation for approx. dynamic programming
%
% Discrete function approximation for use with approximate dynamic
% programming and other simulations. For this class, the state must have
% already been translated (or naturally be) the set of integers between zero and the maximumthe
% range for each dimension
%
% originally by Bryan Palmintier 2011

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
% 0.1  2011-03-17 19:20  BryanP      Initial Version
% 0.2  2011-03-24 12:21  BryanP      Debug multiple state indexing
% 0.3  2011-03-24 13:45  BryanP      It works! (based on manual checking)
% 0.4  2011-03-24 15:45  BryanP      Extracted mat2ind for general use
% 0.5  2011-03-24 20:31  BryanP      Further debugs & new copy() method
% 0.6  2011-03-25 18:21  BryanP      Added raw() method
% 0.7  2012-03-25 15:45  BryanP      Overhauled for FuncApprox v4 
% 0.8  2016-04-25 13:29  HongyuW     To accomodate the new FuncApprox
% 0.9  2017-06-01 10:15  NicolasG    Update bug fix (python backporting)
% 1.0  2017-06-01 14:30  NicolasG    Handle duplicate points in update (python backporting)

%     properties (Access=protected)
%         %Value and meta-data storage as a collection of n-D arrays
%         Value = [];   % current function approximation for discrete item
%         Count = [];   % number of updates for each item
%         Step  = [];   % current step_size for each item
%         Max_per_dim = []; % easy access to the maximum per dimension
%     end
%     
    properties (SetAccess=protected)
        step_size_fn = @ss1overN; %stepsize convergence function
        Step_opt = []; %stepsize function options
        Range = zeros(2,0);  %two row vector of minimum (1,:) and maximum (2,:) extents
        initial_value = 0;
        Step = [];
        Count = [];
        Value = [];
        max_size_per_dim = []
    end
    
    methods
        %======= Constructor =====
        % changed on April 25th 2016
        function obj = faDiscrete(varargin)

             % initialize parser
            p = inputParser;

            % add the function name to the input parser
            p.FunctionName = mfilename;

            % Defaults

            addRequired(p,'max_size_per_dim' ); % TODO fix -> %@(s) assert( isrow(s), 'max_size_per_dim must be a row vector' ) );
            addOptional(p,'pts', []);
            addOptional(p,'vals', []);
            addOptional(p,'step_size_fn', @ss1overN);

            % Parse varargin
            parse( p, varargin{:} )
            pts = p.Results.pts;
            vals = p.Results.vals;

            %Call super class constructor
            obj = obj@FuncApprox(pts, vals);

            obj.max_size_per_dim        = p.Results.max_size_per_dim;
            obj.step_size_fn            = p.Results.step_size_fn;

            obj.NewPts                  = p.Results.pts;
            obj.NewVals                 = p.Results.vals;

            obj.Value = NaN(obj.max_size_per_dim);

            % Initialize value related metadata
            obj.Count = zeros(obj.max_size_per_dim);
            obj.Step = ones(obj.max_size_per_dim);
            % TODO input parser
            % TODO step size fnc
            % TODO initial value
            % TODO add debug to calculate StoreVals

            % TODO Get pts to work with non integers

        end
        
    end
    
    methods (Access = protected)
        
        function build_func(obj, varargin)    
        % This function is used to maintain an internal represention of the
        % lookup table, regression, etc.

            % TODO : Convert non integer states to integer coordinates,
            % which get converted to the index of the Value array

            % TODO : check if it doesn't exist, if it doesn't then expand
            % value, step and count array
            
            % Find linear indicies from the state list

            idx = mat2ind(obj.max_size_per_dim, obj.NewPts);

            % TODO : Check for duplicates in idx
            
            never_visted_mask = (obj.Count(idx) == 0);
            
            %Following code makes sure we do not update the stepsize
            %of points we previously stored but do not observe on the
            %current update (Nicolas-06/01/2017)
            if nargin>=1
               new_points=varargin{1};
               vals=[];
               if nargin>=2
                   vals= varargin{2};
               end
               idx=mat2ind(obj.max_size_per_dim, new_points);
               never_visted_mask = (obj.Count(idx) == 0);
               obj.Count(idx) = 1 + obj.Count(idx);
               obj.Step(idx) = obj.step_size_fn(obj.Count(idx), obj.Step(idx), obj.Step_opt);
            else
                obj.Count(idx) = 1 + obj.Count(idx);
                obj.Step(idx) = obj.step_size_fn(obj.Count(idx), obj.Step(idx), obj.Step_opt);

                vals = obj.NewVals;
            end
           %End of modifications (Nicolas-01/06/2017)
            
           % Put a valid number in all never_visited_mask points to
           % following math will work without NaNs
           obj.Value(idx(never_visted_mask)) = 0;

           % Actually update values weighted by stepsize
           obj.Value(idx) = (1 - obj.Step(idx)) .* obj.Value(idx) + obj.Step(idx).*vals;
        
           %if resulting values, etc. are requested, use approx to find
           %them
           if nargout > 0
                [values, stds, steps] = approx(obj, states);
           end
                        
        end
        
       
        function out_vals = do_approx(obj, out_pts, varargin)             

            % Find linear indicies from the state list
            idx = mat2ind(obj.max_size_per_dim, out_pts);

            display(idx)
            % Extract the values
            out_vals = obj.Value(idx);
        
            %TODO: handle standard deviations (or variances)
            % Standard deviations are not currently stored/computed, so return NaN
            stds = NaN*ones(size(out_vals));
            
            %return current step sizes
            steps = obj.Step(idx);
        
        end
        
    end
    
    methods
        
        function out_vals=update(obj, pts, vals, varargin)
            %Update function for faDiscrete
            %Acts as a wrapper around super class update method
            %The purpose is to handle duplicate points in update queries
            %by splitting the initial query in multiple subqueries without
            %duplicates.
            %See duplicate_row_count.m and update_scheduler.m
            %Nicolas - 06/01/2017
            
            %Call update_scheduler to break in subqueries
            [p_list, v_list]=update_scheduler(pts, vals); 
            
            %Get the number of subqueries
            n_queries=size(p_list,2);
            
            %If that number is equals to one, there is no
            %duplicates. We can call the super class update
            %method the the original arguments
            if n_queries==1
                out_vals=update@FuncApprox(obj, pts, vals);
            %Otherwise, there are duplicates in the original query
            else
                %For each subquery, call the super class update method
                for query_idx=1:n_queries
                    out_vals=update@FuncApprox(obj, p_list{query_idx}, v_list{query_idx});
                end
            end
        end

       
        function test(obj)
           
            obj.build_func()
            
        end
        
    end
    
end

% 
% function [values, stds, steps] = update(obj, states, values)
% % [values, stds, steps] = obj.update(states, values)
% %
% % Upates the value(s) for (a) given state(s). If needs_jacobian is
% % true, callers are expected to pass a valid jacobian. For lists of
% % states, the jacobian is passed as a cell column vector
% 
%     % Find linear indicies from the state list
%     idx = mat2ind(size(obj.Value),states);
% 
%     obj.Count(idx) = 1 + obj.Count(idx);
%     obj.Step(idx) = obj.Step_fun(obj.Count(idx), obj.Step(idx), obj.Step_opt);
%     obj.Value(idx) = (1 - obj.Step(idx)) .* obj.Value(idx) + obj.Step(idx).*values;
% 
%     %if resultilng values, etc. are requested, use approx to find
%     %them
%     if nargout > 0
%         [values, stds, steps] = approx(obj, states);
%     end
% 
% end



