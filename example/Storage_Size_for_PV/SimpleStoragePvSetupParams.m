function params = SimpleStoragePvSetupParams(params_in)
% SIMPLESTORAGEPVPARAMSSETUP Parameter defaults and manipulations for simple storage+PV problems.
%
% params = SimpleStoragePvParamsSetup(params_in)
%    Enables selectively overwriting of defaults by specifying one or more
%    alternate parameter settings as a structure or name-value paired cell
%    array. Required fields:
%         'n_periods'       Number of decision periods--Included b/c we need it for expanding some fields
%         'storage_kVA'     List of valid storage power ratings (kVA) 
%         'storage_kWh'     List of valid storage capacities (kWh)
%
% See also: SimpleStoragePvProblemSetup

% Implementation Notes:
%  - Structures are generally passed by reference, so there is minimal
%    overhead in passing the large multi_inv_prob structure around between
%    functions. for more details see:
%    http://www.mathworks.com/support/solutions/en/data/1-15SO4/index.html?solution=1-15SO4
%
% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   3  2017-07-16 21:24  BryanP      Completed initial parameter, etc. definition 
%   2  2017-07-16 00:31  BryanP      Revised structure of pv_adoption uncertainty for rpTransMatrix 
%   1  2017-07-11 19:27  BryanP      adapted from MultiInvSetupParams v6 with data from FlexDesign.xlsx v13 


%--- Handle defaults (Default to ISGT paper parameters)
default_params = { %Required inputs
                    'n_periods'         []     % Number of decision periods, Included b/c we need it for expanding some fields
                    'storage_kVA'       {}     % List of valid storage power ratings (kVA) over time 
                    'storage_kWh'       {}     % List of valid storage capacities (kWh) over time
                    
                   %System and DER assumptions  
                    'feeder'            'ieee_13'   % String name of feeder to simulate in OpenDSS
                    'days_to_simulate'  {'00-Jan-2000', '00-Apr-2000', '00-Jul-2000', '00-Oct-2000'} %List of days to simulate
                    'run_opendss'       true        % Option to Use OpenDSS for operations simulations. If false uses the stored_ops lookup table.
                    'stored_ops'        'not_used'  % If run_opendss is false, use this lookup table: 3 cols for state, backfeed, 24 hrly net energy 
                    
                    'pv_state_set'  {0.76, [0.76+0.33; 0.76+0.1], [0.76+0.33+0.42; .76+0.1+0.4; 0.76+0.1+0.25]} % Cell vector listing possible states over time in rpTransMatrix format
                    'pv_trans_set'  {[1 1]/2, [0.80 0.20  0   ; 0   0.2 0.8]} % Cell vector listing transition probabilities over time

                    'storage_hr_limits' [1, 5]  % Limit storage power/energy combinations based on hour rating

                   %Operations Assumptions (row vector over time, last entry repeated if needed)
                    'ops_to_annual'     91.5   % Scale factor to scale simulated operations to annual values
                    'ops_backfeed_cost' -270    % Backfeed cost/penalty ($/MWh)
                    'ops_energy_cost'   0       % Wholesale energy supply cost ($/MWh)
                    
                    'ops_pre_decision_yrs'  3   % Number of years to scale pre-decision operations cost by
                    'ops_post_decision_yrs' 1   % Number of years to scale post-decision operations cost by
                    'ops_terminal_yrs'      15  % Number of years to scale terminal period operations cost by
                    
                   %Equipment price assumptions (row vector over time, last entry repeated if needed)
                    % transformer: estimated as power low based on data in Cossi et al., ?Primary Power Distribution Systems Planning Taking into Account Reliability, Operation and Expansion Costs.?
                    %   cost = base * kVA ^ expon
                    'xfmr_cost_base'	-499.4   % $/kVA
                    'xfmr_cost_expon'   0.6638   % unitless
                    
                    'inv_cost'          -[150 90]    % inverter  ($/kVA)
                    'bos_cost'          0           % balance of system ($/kVA)
                    % battery cost (sources: 2015: Arena, 2020: Tesla Powerpack) 
                    'bat_cost'          -[352 250]   % $/KWh
                    
                    'scrap_ratio'               0.2     % Percent of value if scrapped/downsized
                    'replace_on_inv_upsize'     [1 1]   % Map of storage components to replace if inverter increases in size
                    'scrap_value_on_inv_upsize' [1 1]   % Map of storage components with scrap value if inverter increases in size
                 };

params = DefaultOpts(params_in, default_params);

%Confirm required params entered
for fld_idx = 1:size(default_params,1)
    fld_name = default_params{fld_idx, 1};
    if isempty(params.(fld_name))
        error('Adp:MissingRequiredInput', ...
                'The field %s must be speficied (not empty) in either the input parameter structure or as a string value pair', fld_idx)
    end
end

%=== Check sizes: vectorize any lingering scalars and extend/pad if too short 
field_names   = {'storage_kVA', 'storage_kWh', ...
    'ops_backfeed_cost', 'ops_energy_cost', 'xfmr_cost_base', ...
    'xfmr_cost_expon', 'inv_cost', 'bos_cost', 'bat_cost'};

for f = 1:length(field_names)
    % -- Fix bad values
    if isempty(params.(field_names{f})) || (isnumeric(params.(field_names{f})) && any(isnan(params.(field_names{f}))))
        %Set any unreasonable values (back) to the default. Do this
        %first so that scalar defaults will expand properly
        params.(field_names{f}) = default_params.(field_names{f});
        %TODO: add warning
    end

    % -- Pad, duplicate, and trim as needed
    if length(params.(field_names{f})) == 1
        %Expand any scalars to apply to time periods
        % Note: repmat also works with cell vectors
        params.(field_names{f}) = repmat(params.(field_names{f}), 1, params.n_periods + 1);
    else
        %Adjust any existing vectors to be of the correct length:
        if length(params.(field_names{f})) > params.n_periods + 1
            % truncate if too long (also works with cell vectors)
            params.(field_names{f}) = params.(field_names{f})(1:params.n_periods+1);
        elseif length(params.(field_names{f})) < params.n_periods +1
            % extend with final values if too short (also works with cell vectors)
            params.(field_names{f})(length(params.(field_names{f})):(params.n_periods + 1)) = ...
                    params.(field_names{f})(length(params.(field_names{f})));
        end
    end
end
end 