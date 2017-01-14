function multi_inv_problem = MultiInvInit(total_space, prod_space, p_demand, multi_inv_problem_in)
% MULTIINVINIT initialize defaults for multi-product inventory problem (DP & ADP)
%
% Usage: multi_inv_prob = MultiInvInit(total_space, prod_space, p_demand, multi_inv_prob)
%
% Parameters
%      total_space  scalar value of total storage/shelf space shared by all
%                   products (in any unit)
%      prod_space   A row vector of unit space taken by each product
%      p_demand     Probability of demand, as either
%                    + a vector with length(prod_space) elements
%                      interpreted as a poisson lambda value of independant
%                      demands for each product
%                    + a 1-D cell array with an entry for each product
%                      containing a vector of marginal demand probabilities
%                      for the corresponding product. Note that the demand
%                      vectors are assumed to start at demand = 0.
%                    + NOT IMPLEMENTED YET an n-dim array whose elements
%                      sum to one representing the joint pmf of demands for
%                      products. For help on introducing correlation into
%                      such an array see the MATLAB help for copulas by
%                      typing doc copulas
%      mmulti_inv_problem       An optional structure that predefines some or
%                   all of the other inventory functions and parameters to
%                   non-default values.
%
% MultiInventory multi_inv_prob initializes a standard problem structure by
% defining default values for all MultiInv related parameters. These can be
% overridden if desired by passing in an existing multi_inv_prob with these
% already defined.
%  FUNCTIONS (at root level)
%
%  In addition the following functions default to null values, but can be
%  overwritten if needed:
%
%   PROBLEM-SPECIFIC-PARAMETERS (passed in params sub-structure)
%    params:
%      order_cost       fixed cost per order (default=4)
%      unit_cost        per_unit order cost (default=2 for all products)
%      hold_cost        per_unit_per period holding cost (default=1 per unit)
%      sales_price      per unit sales price (default=8 for all products)
%      term_unit_value  value per unit in terminal period (default=0 for all products)
%
% Notes:
%  - this implementation uses extra memory to store the entire list of
%    states, the demand probability matrix, etc. in hopes of faster
%    computation. For large problems, a different representation that uses
%    less storage by dynamically generating states, probabilities, etc. on
%    the fly may be required, with the trade-off of potentially slower
%    computation.
%
% see also:
%   MultiInvCost, MultiInvState, MultiInvTermValue, MultiInvTransProb, MultiInvDecision
%
% originally by Bryan Palmintier 2010
% overhauled for shared problem structure by Bryan Palmintier 2016

% Implementation Notes:
%  - Structures are generallypassed by reference, so there is minimal
%    overhead in passing the large multi_inv_prob structure around between
%    functions. for more details see:
%    http://www.mathworks.com/support/solutions/en/data/1-15SO4/index.html?solution=1-15SO4
%
% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-15 09:47  BryanP      adapted and expanded from InvInit v2
%   2  2010-05-15 16:45  BryanP      added storage of state list
%   3  2010-05-16 01:40  BryanP      added reverse state lookup table
%   4  2010-05-16 21:30  BryanP      converted lookup to n-D array (faster)
%   5  2010-05-26 07:25  BryanP      updated for ValidCombin v4
%   6  2011-03-29 10:10  BryanP      Renamed MultiInvValidCombin to CombinWithLimits
%   7  2016-07-06 23:50  BryanP      Overhauled to use new shared problem structure

%--- Establish defaults

default_params = { %Warehouse and product space, size, & demand configuration
                    'total_space'       NaN     %Input to parent
                    'prod_space'        NaN     %Input to parent
                    'num_products'      NaN     %determined from prod_space
                    'max_inv'           NaN     %determined from prod_space & total_space
                    'demand'            NaN     %a structure built from p_demand
                   %Pricing assumptions
                    'order_cost'        -4
                    'term_unit_val'     0       %will be vectorized as needed
                    'unit_cost'         -2      %will be vectorized as needed
                    'hold_cost'         -1      %will be vectorized as needed
                    'sales_price'       8       %will be vectorized as needed
                 };


multi_inv_defaults = {
                        'setState'            NaN

                        % Returns an set* object given the current pre-state
                        'fDecisionSet'        NaN

                        % A cell vector of RandProc objects
                        'RandomItems'         NaN

                        'faDecisionCost'      NaN

                        'fDecisionApply'      @MultiInvDecisionApply
                        'fUncertainApply'
                        'fTerminalValue'      @MultiInvTermValue

                        'fOpsBeforePre'       @MultiInvSim
                        'fOpsAfterPre'        NaN
                        'fOpsBeforePost'      NaN
                        'fOpsAfterPost'       NaN

                        'fOptimalDecision'    NaN

                        'fSimUncertainty'     NaN
                        'fPostVfunMap'        NaN
                        'fDecCompare'         NaN
                        'fValCompare'         NaN

                        'params'              default_params
                     };

%--- Build up problem structure
%Setup an empty set of inputs if no existing problem passed in.
%That way we can still use SetDefaultOpts
if nargin < 4
    multi_inv_problem_in = {};
end

% Actually merge these external values to create our bare problem strcuture
multi_inv_problem = SetDefaultOpts(multi_inv_problem_in, multi_inv_defaults);

%--- Update problem specific parameters
% Start with any values provided from the existing problem
multi_inv_problem.params = SetDefaultOpts(multi_inv_problem_in.params, default_params);
% Now adjust things for the required parameters
multi_inv_problem.params = ParamsSetupHelper_(total_space, prod_space, p_demand, default_params, multi_inv_problem.params);

%--- Setup Decision Cost Function (Approximation)
decision_intercept = zeros(1, multi_inv_problem.params.num_products);
%Note: unit_cost is set to be a vector in ParamsSetupHelper_()
decision_slope = multi_inv_problem.params.unit_cost;
multi_inv_problem.faDecisionCost = faSepAfineOneTime(decision_intercept, decision_slope);

%--- Additional setup

% TODO: Check existing StateSet rather than rebuilding it


%-- Build up a table (structure) of states and associated parameters --
% Recursively compute a matrix list of all valid multi-product inventory
% states.
%
% Each row corresponds to a state with a column containing the level of
% inventory for each product type
multi_inv_problem.states = CombinWithLimits(multi_inv_problem.prod_space, multi_inv_problem.total_space);

multi_inv_problem.n_states = size(multi_inv_problem.states,1);

%-- Now store a reverse lookup table:
% The goal is to have an easy way to indentify the state index number for a
% given state. This feature is required when computing the transition
% probabilites for a stochastic transition.
%
% Note: may need to uncomment one of these if you play around with
% alternative MultiInvLookup algorithms
% 1) a java hashtable (commented out)
% 2) a MATLAB structure
% 3) just a big n-D array

%multi_inv_prob.state_hash = java.util.Hashtable; % #1 (java hash table) only
% for s = 1:multi_inv_prob.n_states;                % #1 & 2 (MATLAB struct)
%     key= sprintf('i%d',multi_inv_prob.states(s,:));            % #1 & 2
% %    multi_inv_prob.state_hash.put(key, java.lang.Integer(s)); % #2 only
%     multi_inv_prob.state_lookup_struct.(key)=s;                % #1 & 2
% end                                                    % #1 & 2

%#3 only...
% Build up a cell array containing the columns of the state vector to use
% as subscripts into the array.
% Note: add one b/c the 0 inventory level is an invalid MATLAB subscript
subscript_list = cell(1,multi_inv_problem.num_products);
for p = 1:multi_inv_problem.num_products
    subscript_list{p} = multi_inv_problem.states(:,p) + 1;
end

%#3 initialize the n-D array
multi_inv_problem.state_lookup = NaN(multi_inv_problem.max_inv+1);
% And now extract the equivalent single index for all of our states. (#3)
%
% Note: using {:} effectively converts each element of subscript_list
% to a separate parameter, thereby providing one column for each dimension
% as is required by sub2ind.
idx_list = sub2ind(size(multi_inv_problem.state_lookup), subscript_list{:});

%And now populate our n-D array with the state number corresponding to any
%valid state (#3)
multi_inv_problem.state_lookup(idx_list) = 1:multi_inv_problem.n_states;

%-- And also store the possible demands for each product, to avoid having
%to recompute with each call to the transition probability matrix

% first find the max demand with a non-zero probability.
%multi_inv_prob.max_demand = zeros(1,multi_inv_prob.num_products);
for p = 1:multi_inv_problem.num_products;
    multi_inv_problem.demand(p).max = length(multi_inv_problem.demand(p).prob)-1;

    % and store a list of all the possible demand values for each product
    multi_inv_problem.demand(p).list = (0:multi_inv_problem.demand(p).max)';
end

end %main function





function params = ParamsSetupHelper_(total_space, prod_space, p_demand, default_params, params_in)
% Helper function to manage defaults and othe manipulations for
% Multi-Inventory problems.
    % --- Now overide any existing values with those explicitly required:
    params.total_space = total_space;

    %Store the space requirements of each product for future use
    params.prod_space = prod_space;

    % Extract the number of products from the product space vector
    params.num_products = length(prod_space);

    %Determine the maximum inventory of each product, if it were to occupy the
    %entire shelf space, this allows us to limit the demand space for poissons,
    %and may get used elsewhere, too.
    params.max_inv = floor(total_space ./ prod_space);

    % --- Handle Demand Probabilities
    if iscell(p_demand)
        % Cell contains specified vectors of demand per product
        %TODO add dimension checks
        for p = 1:params.num_products;
            %ensure we have a column for consistancy with the poisspdf
            %representation for vectorized operations
            if isrow(p_demand{p})
                p_demand{p} = p_demand{p}';
            end
            params.demand(p).prob = p_demand{p};
        end
    else
        [rows, cols] = size(p_demand);
        if cols == 1 || rows == 1
            % Treat vectors as a set of lambdas for independant poissons
            assert(length(p_demand) == params.num_products, 'Must specify a poisson lambda demand for each item')
            for p = 1:params.num_products;
                params.demand(p).prob = poisspdf(0:(params.max_inv(p)-1),p_demand(p))';
                % lump any additional demand into the
                % highest demand... this seems consistant with the idea that any
                % additional demand can only produce the same ammount of sales as the
                % maximum inventory
                params.demand(p).prob(params.max_inv(p)+1) = 1 - sum(params.demand(p).prob);
            end
        elseif (ndims(p_demand) == params.num_products && params.num_products >=2) ...
                || (ndims(p_demand) ==1 && params.num_products == 1)
            % If we have the right number of dimensions in p_demand treat as joint pmf
            error('Support for correlated demand probabilities not yet implemented')
            %multi_inv_prob.p_demand = p_demand;
        else
            error('Invalid demand array.')
        end
    end

    %check sizes and vectorize any lingering scalars
    field_names   = {'term_unit_val', 'unit_cost', 'hold_cost','sales_price'};

    for f = 1:length(field_names)
        if isempty(params.(field_names{f})) || isnan(params.(field_names{f}))
            %Set any unreasonable values (back) to the default. Do this
            %first so that scalar defaults will expand properly
            params.(field_names{f}) = default_params.(field_names{f});
        end

        if length(params.(field_names{f})) == 1
            %Expand any scalars to apply to all products
            params.(field_names{f}) = ones(1,params.num_products)*params.(field_names{f});
        else
            %Adjust any existing vectors to be of the correct length:
            if length(params.(field_names{f})) > params.num_products
                % truncate if too long
                params.(field_names{f}) = params.(field_names{f})(1:length(params.(field_names{f})));
            elseif length(params.(field_names{f})) < params.num_products
                % pad with defaults if too short
                params.(field_names{f})(length(params.(field_names{f})):params.num_products) = ...
                        default_params.(field_names{f});
            end
        end
    end
end  %ParamsSetupHelper_
