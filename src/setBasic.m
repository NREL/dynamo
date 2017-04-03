classdef setBasic < AbstractSet
%setBASIC generic set for mix of discrete and continuous rectalinear spaces 
%
% see also:
%  setCombinWithLimits, IntegerRangeFromReal, RealFromIntegerRange 
%
% adapted from SampleNdRange (v4) by Bryan Palmintier 2016

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   6  2016-10-21 15:24  BryanP      BUGFIX: correct handling of multi-dimensional ranges 
%   5  2016-10-06 11:24  BryanP      Fix bug by explicitly calling superclass constructor 
%   4  2016-09-29 14:35  BryanP      Fix bug with repeating creation of qrand stream... oops 
%   3  2016-09-29 14:10  BryanP      Updated to match AbstractSet v2 (and v3) 
%   2  2016-07-08 02:20  BryanP      Various bug fixes. Now working 
%   1  2016-07-07 22:20  BryanP      Adapted from SampleNdRange v4 with selected pieces from RandProc v5 

    properties
        DStepsForContinuous = 10;  %Number of uniform divisions for continuous dimenstion with dlist
    end

    % Read-only properties (set by constructor)
    properties (SetAccess='protected')
        MinMax                  % Range per dimension as 2 rows [min; max]
        StepSize
    end
        
    % Hidden properties
    properties (Access='protected')
        FullValueList = [];
    end
    
    methods
        function obj = setBasic(sample_minmax, sample_type, sample_step)
        % CONSTRUCTOR (setBasic)
        %
        % obj = setBasic(sample_minmax, sample_type)
        %
        % obj = setBasic(sample_minmax, sample_type, sample_step)

            %Support blank calls per object requirements for serializing
            if nargin == 0
                sample_minmax = [];
            end

            if nargin < 2 || isempty(sample_type)
                sample_type = 'rand';
            end
            
            %IMPORTANT: must call superclass to match proper inputs
            % Also sets N_dim
            obj@AbstractSet(size(sample_minmax,2), sample_type);

            if nargin < 3 || isempty(sample_step)
                sample_step = zeros(1, obj.N_dim);
            end
            
            
            obj.MinMax = sample_minmax;
            %Compute the sample range between min (first row) and max
            %(second row) using range = max-min
            obj.SampleParams.Range = diff(obj.MinMax,1, 1);
            obj.SampleParams.Offset = obj.MinMax(1,:);
            
            %If scalar step, assume it applys to all dimensions
            if isscalar(sample_step)
                sample_step = ones(1, obj.N_dim) * sample_step;
            end
            obj.StepSize = sample_step;
            %Set up a boolean mask for any discrete variables
            obj.DiscreteMask = (obj.StepSize ~= 0);

            if any(obj.DiscreteMask)
                %Compute rounded space for discrete values
                [obj.SampleParams.Range(obj.DiscreteMask), obj.SampleParams.Offset(obj.DiscreteMask)] = ...
                    IntegerRangeFromReal(sample_minmax(:, obj.DiscreteMask), obj.StepSize(obj.DiscreteMask));
            end

        end

        function value_list = sample(obj, N, sample_type, varargin) 
        % SAMPLE Return the requested number of multi-dimensional samples
        %  Values are scaled to fill the entire range. Assumed constant
        %  over all time
        %
        % out = setBasicObject.sample()  
        %           Return a single sample
        % out = setBasicObject.sample(N)  
        %           Return N samples
        % out = setBasicObject.sample(N, sample_type)  
        %           First initializes sampling to the specified type (e.g.
        %           'sobol') then returns N samples.

            if nargin < 2 || isempty(N)
                N = 1;
            end
            
            % Attempt to initialize sample if needed
            if nargin >= 3 && not(isempty(sample_type))
                obj.sample_init(sample_type)
            end
            
            % get required samples over 0:1 range using base class
            value_list = obj.SampleParams.fSample(N);
            
            % -- Transform from [0,1] to desired sample_minmax
            % Expand to cover desired rang
            value_list = bsxfun(@times, value_list, obj.SampleParams.Range);
            
            % Round any discrete values
            value_list(:,obj.DiscreteMask) = bsxfun(@times, ...
                floor(value_list(:,obj.DiscreteMask)), obj.StepSize(:,obj.DiscreteMask));
            % Adjust sample upward
            value_list = bsxfun(@plus, value_list, obj.SampleParams.Offset);

        end

        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.
        
        function value_list = as_array (obj, discrete_steps)
        % AS_ARRAY List possible (discrete) states
        %
        % value_list = setBasicObject.as_array()
        %    List all possible states. For continuous values, discretized
        %    using the default number of steps (normally 10)
        %
        % value_list = setBasicObject.as_array(discrete_steps)
        %    Use the specified number of steps, for continuous
        %    approximation

            % If we haven't done so before, build up the full value list
            % (Potentially Slow)
            if isempty(obj.FullValueList)
                % First compute integer range, etc. for all dimensions
                [integer_max, offset, step_size] ...
                    = IntegerRangeFromReal( obj.MinMax, obj.StepSize, obj.DStepsForContinuous );
                
                % Now get list of values for each dimension
                [~, real_values] = RealRangeFromInteger(integer_max, offset, step_size);
                
                % And finally compute all combinations
                obj.FullValueList = allcomb(real_values{:});
            end
            
            %Once built, simply return cached value list
            value_list = obj.FullValueList;
        end
            

        %% ===== General (discrete or continuous) Methods
        
        % RANGE Find value range for given time
        %
        % Returns array with [min; max] value range (assumed constant over
        % time
        function [value_range] = range(obj, ~)
            value_range = obj.MinMax;
        end
            
    end %public methods
    
end

