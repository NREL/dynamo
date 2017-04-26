classdef AbstractSet < handle
%ABSTRACTSET Abstract class for managing samplable multi-timeperiod-aware sets for adp toolbox 
%
% Defines an abstract class (defines the strucuture for related
% subclasses) for a Set use in dynamic programming and
% other simulations
%
% adapted from RandProc by Bryan Palmintier 2016

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   4  2016-09-30 11:13  BryanP      Adding optional pt_dim_names property 
%   3  2016-09-29 14:30  BryanP      Bug fixes from testing with setBasic and demoSets 
%   2  2016-09-29 11:50  BryanP      Updated for only one time period and to allow on the fly sample_setup 
%   1  2016-07-07 09:50  BryanP      Merge Outline from SampleNdRange v4 with selected pieces from RandProc v5 

    properties
        name = '';   %optional name
        pt_dim_names = {};    %optional cell array of names (one item per dimension)       
    end
    
    % Read-only properties (set by constructor)
    properties (SetAccess='protected')
        SampleType = '';    % Must support 'rand' & 'sobol' and gracefully handle unrecognized types
        N_dim = 0;                  % Number of dimensions in set's space
        DiscreteMask;           % Boolean mask of descrete dimesions of set
    end
    
    % Internal properties
    properties (Access='protected')
        SampleParams = struct(  'fSample',     [] ... % (Quasi-)random sample function
                             );  %Note Set and Offset added when using Sobol samples
    end
    
    methods
        %% ===== Constructor
        function obj = AbstractSet(n_dim, sample_type)
        % CONSTRUCTOR Example
        %
        % Abstract example, simply demonstrates sample method initalization
        % sample_type should be supported by all Sets, but n_dim need not
        % be passed at contruction time

            if nargin == 0  %Support no argument constructor syntax using defaults
                sample_type = 'rand';
                n_dim = 1;
            end
            if nargin < 2 || isempty(sample_type)
                sample_type = 'rand';
            end

            %Typically derived from other input in subclasses
            obj.N_dim = n_dim;
            
            obj.sample_init(sample_type);
        end
        
        %% ===== General Sampling Support
        % SAMPLE_INIT_ Initalizes data for basic sampling types
        function sample_init(obj, sample_type, varargin)
            % TODO: support a force re-init
            if strcmpi(obj.SampleType, sample_type)
                %Skip if already initialized
                warning('AdpSample:AlreadyInitialized', 'Sampling already intialized for type "%s". Ignoring', sample_type)
                return
            end
            
            %Otherwise setup sampling as needed
            obj.SampleType = sample_type;

            switch lower(sample_type)
                case 'rand'
                    obj.SampleParams.fSample = @(N) rand(N, obj.N_dim);
                case 'sobol'
                    obj.SampleParams.Set = sobolset(obj.N_dim);
                    obj.SampleParams.Set.scramble('MatousekAffineOwen');

                    % Setup function to provide a series of quasi-random state samples
                    my_stream = qrandstream(obj.SampleParams.Set);
                    obj.SampleParams.fSample = @(N) qrand(my_stream, N);
                otherwise
                    warning('AdpSample:TypeUnknown', 'Unknown sample type %s, defaulting to rand',sample_type)

                    % Default to standard pseudorandom samples & define
                    % with recursion
                    sample_init_(obj, 'rand', varargin) 
                    
            end
        end
        
        
        function value_list = sample(obj, N, sample_type, varargin) 
        % SAMPLE Return the requested number of multi-dimensional samples
        %  AbstractSet implementation simply returns items over a 0:1 range
        %  for each dimension, demonstrates sobol 
        %  sampling and other semantics
        %
        % out = SetObject.sample()  
        %           Return a single sample
        % out = SetObject.sample(N)  
        %           Return N samples using the current (last) sampling type
        % out = SetObject.sample(N, sample_type)  
        %           First initializes sampling to the specified type (e.g.
        %           'sobol') then returns N samples.
            if nargin < 2 || isempty(N)
                N = 1;
            end
            
            % Attempt to initialize sample if needed
            if nargin >= 3 && not(isempty(sample_type))
                obj.sample_init(sample_type)
            end
            
            % simply returns sample over 0:1 (Note "hard part" is in
            % sample_init_) 
            value_list = obj.SampleParams.fSample(N);

        end
        
    end % default methods
    
    methods (Abstract)
        
        %% ===== Support for discrete usage
        % These need to be defined even for continuous processes, for
        % compatability with DP.
  
        % AS_ARRAY List possible (discrete) states
        %
        % value_list = SetObject.as_array()
        %    List all possible states. For continuous values, discretized
        %    using the default number of steps (normally 10)
        %
        % value_list = SetObject.as_array(discrete_steps)
        %    Use the specified number of steps, for continuous
        %    approximation
        value_list = as_array (obj, discrete_steps, varargin)
        

        
        %% ===== General (discrete or continuous) Methods
        
        % RANGE Find value range for given time
        %
        % Returns vector with [min max] value range for specified time
        % if t is not provided, the range for the current simulation time
        % is returned.
        [value_range, state_n_range] = range(obj, t)            
    end %Abstract methods

end

