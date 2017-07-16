classdef setList < AbstractSet
%setLIST Simple set with pre-specified list of options 
%
% The list of members is defined with one state per row in a column vector
% 
% Examples:
%
% >> rng('default'); my_set = setList([1 2; 11 12; 21 22]);
% >> my_set.N_dim
%
% ans =
%   
%        2
%
% >> my_set.sample()
% 
% ans =
%     21    22
% 
% >> my_set.sample(4)
% 
% ans =
%     21    22
%      1     2
%     21    22
%     11    12
% 
% >> my_set.as_array()
% 
% ans =
%      1     2
%     11    12
%     21    22
% 
% 
% see also:
%  setBasic, setCombinWithLimits, AbstractSet 
%
% adapted from setBasic (v6) by Bryan Palmintier 2017

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-07-16 10:15  BryanP      Adapted from setSingleItem v2 

    % Read-only properties (set by constructor)
    properties (SetAccess='protected')
        Values = {};    % one state per row in a column vector
    end
    
    methods
        function obj = setList(values, ~)
        % CONSTRUCTOR (setSingleItem)
        %
        % obj = setNBasic(sample_minmax, sample_type)
        %
        % obj = setBasic(sample_minmax, sample_type, sample_step)

            %Support blank calls per object requirements for serializing
            if nargin == 0
                values = [];
            end
            
            %IMPORTANT: must call superclass to match proper inputs
            % Also sets N_dim
            obj@AbstractSet(size(values, 2));          
            
            obj.Values = values;
            
            %Force only a single dimensional sample, since we just pick
            %from our list
            obj.n_sample_dims = 1;
        end

        function value_list = sample(obj, N, varargin) 
        % SAMPLE Return the requested number of multi-dimensional samples
        %
        % out = setSingleItemObject.sample()  
        %           Randomly sample a single multi-d state
        % out = setSingleItemObject.sample(N)  
        %           Return N samples

            if nargin < 2 || isempty(N)
                N = 1;
            end
            
            rand_idx = randi(size(obj.Values, 1), N, 1);
            value_list = obj.Values(rand_idx, :);
        end

        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.
        
        function value_list = as_array (obj, varargin)
        % AS_ARRAY List possible states (here only one)
        %
        % value_list = setSingleItemObject.as_array()
        %    List all possible states 
        
            %Once built, simply return cached Values list
            value_list = obj.Values;
        end
            

        %% ===== General (discrete or continuous) Methods
        
        % RANGE Find Values range
        %
        % Returns array with [min; max] Values range (assumed constant over
        % time. 
        function [value_range] = range(obj, ~)
            value_range = [ min(obj.Values)
                            max(obj.Values) ];
        end
            
    end %public methods
    
end

