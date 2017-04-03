classdef setSingleItem < AbstractSet
%setSINGLEITEM degenerate set with a single, constant item 
%
% Examples:
%
% >> my_set = setSingleItem([1 2 3 4 5.6]);
% >> my_set.N_dim
%
% ans =
%   
%        5
%
% >> my_set.sample()
% 
% ans =
% 
%     1.0000    2.0000    3.0000    4.0000    5.6000
% 
% >> my_set.sample(4)
% 
% ans =
% 
%     1.0000    2.0000    3.0000    4.0000    5.6000
%     1.0000    2.0000    3.0000    4.0000    5.6000
%     1.0000    2.0000    3.0000    4.0000    5.6000
%     1.0000    2.0000    3.0000    4.0000    5.6000
% 
% >> my_set.as_array()
% 
% ans =
% 
%     1.0000    2.0000    3.0000    4.0000    5.6000
% 
%
% %Note: supports automatic column to row conversions
% >> my_set = setSingleItem([10; 20; 30]);
% >> my_set.Value
% 
% ans =
% 
%     10    20    30
%
%
% see also:
%  setBasic, setCombinWithLimits, AbstractSet 
%
% adapted from setBasic (v6) by Bryan Palmintier 2017

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   2  2017-04-03 10:00  BryanP      Added doctests 
%   1  2017-04-03 09:30  BryanP      Adapted from setBasic v6 

    % Read-only properties (set by constructor)
    properties (SetAccess='protected')
        Value
    end
    
    methods
        function obj = setSingleItem(value_row_vector, ~)
        % CONSTRUCTOR (setBasic)
        %
        % obj = setNBasic(sample_minmax, sample_type)
        %
        % obj = setBasic(sample_minmax, sample_type, sample_step)

            %Support blank calls per object requirements for serializing
            if nargin == 0
                value_row_vector = [];
            end
            
            if iscolumn(value_row_vector)
                value_row_vector = value_row_vector';
            end

            %IMPORTANT: must call superclass to match proper inputs
            % Also sets N_dim
            obj@AbstractSet(length(value_row_vector));          
            
            obj.Value = value_row_vector;
        end

        % Replace some underlying function with dummy shells
        function sample_init(obj, varargin) %#ok<INUSD>
            return
        end
        
        
        function value_list = sample(obj, N, varargin) 
        % SAMPLE Return the requested number of multi-dimensional samples
        %  Here this just duplicates our single value
        %
        % out = setSingleItemObject.sample()  
        %           Returns our value
        % out = setSingleItemObject.sample(N)  
        %           Return N identical samples

            if nargin < 2 || isempty(N)
                N = 1;
            end
            
            value_list = repmat(obj.Value, N, 1);
        end

        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.
        
        function value_list = as_array (obj, varargin)
        % AS_ARRAY List possible states (here only one)
        %
        % value_list = setSingleItemObject.as_array()
        %    List all possible states by returning our single state 
        
            %Once built, simply return cached value list
            value_list = obj.Value;
        end
            

        %% ===== General (discrete or continuous) Methods
        
        % RANGE Find value range for given time
        %
        % Returns array with [min; max] value range (assumed constant over
        % time. Here min=max b/c single value
        function [value_range] = range(obj, ~)
            value_range = repmat(obj.Value, 2, 1);
        end
            
    end %public methods
    
end

