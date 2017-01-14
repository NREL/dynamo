classdef setCombinWithLimits < AbstractSet
%setCOMBINGWITHLIMITS set class for volume restricted Nd spaces
%
% Provides set standard support for listing, sampling, and range
% for volume restricted Nd spaces such as state spaces for multi-product
% inventory or generation capacity expansion planning.
%
% adapted from SampleNdRange (v4) by Bryan Palmintier 2016
%
% see also:
%  setNdRange, AbstactSet, CombinWithLimits

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2016-07-08 07:00  BryanP      Initial version adapted from setNdRange v2 

    properties
        UseValueSpace = true;       %If false returns integer counts
    end

    % Read-only properties (set by constructor)
    properties (SetAccess='protected')
        Opt = struct()      %Structure with all parameters used by CombinWithLimits function
    end
        
    % Hidden properties
    properties (Access='protected')
        FullValueList;
        FullListIsStoredAsValue;
    end
    
    methods
        function obj = setCombinWithLimits(sample_type, use_value_space, varargin)
        % CONSTRUCTOR (setCombinWithLimits)
        %
        % Note: after sample_type, all parameters follow the same
        % conventions as CombinWithLimits. (In fact they are directly passed
        % there)
        %
        % obj = setCombinWithLimits(sample_type, value_space, space_per_item, total_space)
        %   Specify total volume with other defaults
        %
        % obj = setCombinWithLimits(sample_type, value_space, space_per_item, [], item_max)
        %   Specify limits per item (can also be used with total_space)
        %
        % obj = setCombinWithLimits(sample_type, value_space, space_per_item, ...
        %                       total_space, item_max, min_space_used, ...
        %                       item_min, fRoundMax, fRoundMin)

            if nargin < 3  %Support no argument constructor syntax using defaults
                sample_type = 'rand';
                use_value_space = true;
                varargin = {};
            end

            obj.UseValueSpace = use_value_space;    %Note: this can be changed later by user
            obj.FullListIsStoredAsValue = use_value_space;
            
            if obj.UseValueSpace
                [~, obj.FullValueList, ~, obj.Opt] = CombinWithLimits(varargin);
            else
                [obj.FullValueList, ~, ~, obj.Opt] = CombinWithLimits(varargin);
            end
            
            if strcmp(sample_type, 'sobol')
                warning('AdpSample:Unsupported', 'sobol is currently not supported for setCombinWithLimits, defaulting to rand')
                sample_type = 'rand';
            end
            obj.sample_init_(sample_type);
        end

        function value_list = sample(obj, ~, N, varargin) 
        % SAMPLE Return the requested number of multi-dimensional samples
        %  Values are scaled to fill the entire range. Assumed constant
        %  over all time
        %
        % out = NdSampleObject.sample(t)  
        %           Return a single sample
        % out = NdSampleObject.sample(t, N)  
        %           Return N samples
            if nargin < 3 || isempty(N)
                N = 1;
            end
            
            % get required samples over 0:1 range using base class
            value_list = obj.SampleParams.fRaw(N);
            
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
        %
        % IMPORTANT: all processes must have a single, starting state for t=0
        
        function value_list = dlist (obj, ~)
        % DLIST List possible discrete states
        %
        % List possible discrete states (same for all time)

            %If our current Value state is the same as stored, return it
            if obj.UseValueState == obj.FullListIsStoredAsValue
                value_list = obj.FullValueList;
            elseif obj.UseValueState %must be stored as counts
                value_list = bsxfun(@times, obj.FullValueList, obj.Opt.space_per_item);
            else %must need counts from values
                value_list = round(bsxfun(@rdivide, obj.FullValueList, obj.Opt.space_per_item));
            end
            
        end
            

        %% ===== General (discrete or continuous) Methods
        
        % RANGE Find value range (assumed constant over time)
        %
        % Returns array with [min; max] value range (assumed constant over
        % time
        function [value_range] = range(obj, ~)
            %Start with the item count range
            value_range = [ obj.Opt.item_min
                            obj.Opt.item_max ];
                        
            if obj.UseValueSpace
                %and convert to values if needed
                value_range = bsxfun(@times, value_range, opj.Opt.space_per_item);
            end
        end
            
    end %public methods
    
end

