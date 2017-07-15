classdef setCombinWithLimits < AbstractSet
%setCOMBINGWITHLIMITS set class for volume restricted Nd spaces
%
% Provides set standard support for listing, sampling, and range
% for volume restricted Nd spaces such as state spaces for multi-product
% inventory or generation capacity expansion planning.
%
% The set can be represented as either item counts (such as for inventory
% problems) or using the "space" per item, which allows capacity based
% entries such as the installed capacity of power generation.
%
% Examples:
% % Item count space
% >> rng('default')  %Reset random number stream
% >> my_set = setCombinWithLimits('from_list', false, [2 3], 6)
% 
% my_set = 
% 
%   setCombinWithLimits with properties:
% 
%     UseValueSpace: 0
%               Opt: [1×1 struct]
%          N_Combin: 7
%              name: ''
%      pt_dim_names: {}
%        SampleType: 'from_list'
%             N_dim: 2
%      DiscreteMask: []
% 
% >> my_set.as_array()
% 
% ans =
% 
%      0     0
%      0     1
%      0     2
%      1     0
%      1     1
%      2     0
%      3     0
% 
% >> my_set.sample()
% 
% ans =
% 
%      2     0
% 
% >> my_set.sample(5)
% 
% ans =
% 
%      3     0
%      0     0
%      3     0
%      1     1
%      0     0
% 
% 
% % Item count space
% >> my_set = setCombinWithLimits('from_list', true, [2 3], 6)
% 
% my_set = 
% 
%   setCombinWithLimits with properties:
% 
%     UseValueSpace: 1
%               Opt: [1×1 struct]
%          N_Combin: 7
%              name: ''
%      pt_dim_names: {}
%        SampleType: 'from_list'
%             N_dim: 2
%      DiscreteMask: []
% 
% >> my_set.as_array()
% 
% ans =
% 
%      0     0
%      0     3
%      0     6
%      2     0
%      2     3
%      4     0
%      6     0
%
%
% Initially by Bryan Palmintier 2016
%
% see also:
%  setBasic, AbstactSet, CombinWithLimits

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   3  2017-07-14 20:02  BryanP      Update doctests for new output style 
%   2  2017-04-26 04:46  BryanP      Overhauled for revised AbstractSet (v4) 
%   1  2016-07-08 07:00  BryanP      Initial version adapted from setNdRange v2 

    properties
        UseValueSpace = true;       %If false returns integer counts
    end

    % Read-only properties (set by constructor)
    properties (SetAccess='protected')
        Opt = struct();      %Structure with all parameters used by CombinWithLimits function
        N_Combin = NaN;
    end
        
    % Hidden properties
    properties (Access='protected')
        FullValueList;
        FullListIsStoredAsValue;
    end
    
    methods
        function obj = setCombinWithLimits(sample_type, use_value_space, space_per_item, varargin)
        % CONSTRUCTOR (setCombinWithLimits)
        %
        % obj = setCombinWithLimits(sample_type, use_value_space, space_per_item, total_space)
        %   Specify the sampling type: 'rand' or 'sobol' for accept/reject [NOT IMPLEMENTED] 
        %                              'from_list' to build full set of combin and sample from the list 
        %   use_value_space (boolean): return combinations as total 
        %                              capacity/space (true) or item counts (false)
        %   remaining options follow the same conventions as
        %   CombinWithLimits. (passed directly). Here is the most common
        %   use with a finite space_per_item used to fill a total_space
        %
        % obj = setCombinWithLimits(sample_type, value_space, space_per_item, [], item_max)
        %   Specify limits per item (can also be used with total_space)
        %
        % obj = setCombinWithLimits(sample_type, value_space, space_per_item, ...
        %                       total_space, item_max, min_space_used, ...
        %                       item_min, fRoundMax, fRoundMin)
        %   Full CombinWithLimits options. See help there

            if nargin < 3  %Support no argument constructor syntax using defaults
                sample_type = 'rand';
                use_value_space = true;
                varargin = {};
            end

            if strcmp(sample_type, 'from_list')
                %Explicitly build combin with limits list and sample from
                %it. Should be fastest for higher dimensions and thin
                %shells
                %
                % For now this is the default
                obj.SampleType = sample_type;
            else
                %TODO: Setup accept/reject sampling using sobol, multi-dimensional rand 
                % obj@AbstractSet(size(space_per_item, 2), sample_type);
                if strcmp(sample_type, 'sobol')
                    warning('AdpSample:Unsupported', 'sobol is currently not supported for setCombinWithLimits, defaulting to from_list')
                end
                %For rand, use from_list without warning
                obj.SampleType = 'from_list';
            end
            
            %Update N_dim b/c we don't call parent constructor
            obj.N_dim = size(space_per_item, 2);
            obj.UseValueSpace = use_value_space;    %Note: this can be changed later by user
            obj.FullListIsStoredAsValue = use_value_space;
            
            % If needed, build the full list using CombinWithLimits
            if obj.UseValueSpace
                % The second return value includes the capacity-based list
                [~, obj.FullValueList, ~, obj.Opt] = CombinWithLimits(space_per_item, varargin{:});
            else
                % The first return value includes item counts
                [obj.FullValueList, ~, ~, obj.Opt] = CombinWithLimits(space_per_item, varargin{:});
            end
            obj.N_Combin = size(obj.FullValueList, 1);
            
        end

        function value_list = sample(obj, N, ~, varargin) 
        % SAMPLE Return the requested number of multi-dimensional samples
        %  Values are scaled to fill the entire range. Assumed constant
        %  over all time
        %
        % out = setCombinWithLimitsObject.sample()  
        %           Return a single sample
        % out = setCombinWithLimitsObject.sample(N)  
        %           Return N samples
            if nargin < 2 || isempty(N)
                N = 1;
            end
            
            if not(strcmp(obj.SampleType, 'from_list'))
                warning('AdpSample:Unsupported', 'Unsupported sample type (%s) is currently not supported for setCombinWithLimits, defaulting to from_list', obj.SampleType)
            end
            
            %TODO: implement sample-reject
            
            sample_idx = randi(obj.N_Combin, N, 1);
            value_list = obj.FullValueList(sample_idx, :);
        end

        %% ===== Support for discrete usage
        
        function value_list = as_array (obj, ~, varargin)
        % AS_ARRAY List possible (discrete) states
        %
        % setCombinWithLimit is inherently already discrete, so simply
        % return our list of item combinations

            %If our current Value state is the same as stored, return it
            value_list = obj.FullValueList;
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

