function params = MultiInvParamsSetup(params_in)
% MULTIINVPARAMSSETUP Manage parameter defaults and othe manipulations for
% Multi-Inventory problems.
%
% params = MultiInvParamsSetup()
%    Setups a complete multi-inventory problem-specific parameter
%    structure using all default settings. In addition to raw parameter
%    settings, scalar values are converted to vectors as needed, and
%    probabilities are autodefined if not specified.
%
% params = MultiInvParamsSetup(params_in)
%    Enables selectively overwriting of defaults by specifying one or more
%    alternate parameter settings as a structure or name-value paired cell
%    array.
%
% 
% See also: MultiInvProblemSetup

% Implementation Notes:
%  - Structures are generally passed by reference, so there is minimal
%    overhead in passing the large multi_inv_prob structure around between
%    functions. for more details see:
%    http://www.mathworks.com/support/solutions/en/data/1-15SO4/index.html?solution=1-15SO4
%
% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   6  2017-06-14 05:08  BryanP      Expanded documentation 
%   5  2017-04-26 17:32  BryanP      Add support for single item 
%   4  2016-10-21 16:55  BryanP      Bug fix: handle vector values in checks for cost parameters 
%   3  2016-10-07 01:24  BryanP      Bug fix: actually pass poisson probability vectors to p_demand 
%   2  2016-10-06 11:24  BryanP      Finished clean-up. Seems to work 
%   1  2016-10-02        BryanP      extracted and adapted from MultiInvInit v7 


%--- Handle defaults
default_params = { %Warehouse and product space, size, & demand configuration
                   %Required inputs
                    'total_space'       NaN     %Max space in warehouse
                    'unit_space'        NaN     %Space per item
                    'p_demand'          NaN     %Probability per item, poisson if scalar per item, cell with columns otherwise
                   %Pricing assumptions
                    'order_cost'        -4
                    'term_unit_val'     0       %will be vectorized as needed
                    'unit_cost'         -2      %will be vectorized as needed
                    'hold_cost'         -1      %will be vectorized as needed
                    'sales_price'       8       %will be vectorized as needed
                   %Computed here (so not included yet, so can check fields)
                   % 'n_products'      NaN     %determined from unit_space
                   % 'max_inv'           NaN     %determined from unit_space & total_space
                 };

params = DefaultOpts(params_in, default_params);

%Confirm required params entered
for fld_idx = 1:size(default_params,1)
    fld_name = default_params{fld_idx, 1};
    if isnan(params.(fld_name))
        error('Adp:MissingRequiredInput', ...
                'The field %s must be speficied (not NaN) in either the input parameter structure or as a string value pair', fld_idx)
    end
end

%=== Compute additional fields
% Extract the number of products from the product space vector
params.n_products = length(params.unit_space);

%Determine the maximum inventory of each product, if it were to occupy the
%entire shelf space, this allows us to limit the demand space for poissons,
%and may get used elsewhere, too.
params.max_inv = floor(params.total_space ./ params.unit_space);

% --- Handle Demand Probabilities
if iscell(params.p_demand)
    % Cell contains specified vectors of demand per product
    %TODO add dimension checks
    for p = 1:params.n_products;
        %ensure we have a column for consistancy with the poisspdf
        %representation for vectorized operations
        if isrow(params.p_demand{p})
            params.p_demand{p} = params.p_demand{p}';
        end
    end
else
    [rows, cols] = size(params.p_demand);
    if  params.n_products == 1 && length(params.p_demand) > 1
        %Single item with specified probabilities
        if length(params.p_demand) < params.max_inv + 1
            %Pad an unspecified demand with zero
            params.p_demand(params.max_inv + 1) = 0;
        end
        params.p_demand = {params.p_demand};
    elseif cols == 1 || rows == 1
        % Treat vectors as a set of lambdas for independant poissons
        assert(length(params.p_demand) == params.n_products, 'Must specify a poisson lambda demand for each item')
        temp_p_demand = cell(1, params.n_products);
        
        for p = 1:params.n_products;
            temp_p_demand{p} = poisspdf(0:(params.max_inv(p)-1),params.p_demand(p))';
            % lump any additional demand into the
            % highest demand... this seems consistant with the idea that any
            % additional demand can only produce the same ammount of sales as the
            % maximum inventory
            temp_p_demand{p}(params.max_inv(p)+1) = 1 - sum(temp_p_demand{p});
        end
        params.p_demand = temp_p_demand;
    elseif (ndims(params.p_demand) == params.n_products && params.n_products >=2) ...
            || (ndims(params.p_demand) ==1 && params.n_products == 1)
        % If we have the right number of dimensions in p_demand treat as joint pmf
        error('Support for correlated demand probabilities not yet implemented')
    else
        error('Invalid demand array.')
    end
end

    %check sizes and vectorize any lingering scalars
    field_names   = {'term_unit_val', 'unit_cost', 'hold_cost','sales_price'};

    for f = 1:length(field_names)
        if isempty(params.(field_names{f})) || any(isnan(params.(field_names{f})))
            %Set any unreasonable values (back) to the default. Do this
            %first so that scalar defaults will expand properly
            params.(field_names{f}) = default_params.(field_names{f});
            %TODO: add warning
        end

        if length(params.(field_names{f})) == 1
            %Expand any scalars to apply to all products
            params.(field_names{f}) = ones(1,params.n_products)*params.(field_names{f});
        else
            %Adjust any existing vectors to be of the correct length:
            if length(params.(field_names{f})) > params.n_products
                % truncate if too long
                params.(field_names{f}) = params.(field_names{f})(1:length(params.(field_names{f})));
            elseif length(params.(field_names{f})) < params.n_products
                % pad with defaults if too short
                params.(field_names{f})(length(params.(field_names{f})):params.n_products) = ...
                        default_params.(field_names{f});
            end
        end
    end
end 