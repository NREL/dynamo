function  params  = HEMS_setup_params( params_in )
%function [ output_args ] = HEMS_setup_params( input_args )
%   This functions sets up the HEMS parameters to as is used in the DYNAMO
%   toolkit. 
default_params = {
                % See HEMS_Main for the explanation of these parameters.
                'time_slots'            4
                'total_space'           realmax('single') 
                'unit_space'            [1] 
                'appliance_range'       {[70 80]} 
                'stochastic_range'      {[75,85]} 
                'stochastic_medi'       []  
                'stochastic_sigma'        [0.1]
                'stochastic_prob'       0    %calculate later.
                'unit_cost'            [-2] 
                'desired_value'        [75] 
                'arrive_departure_time' [1,4] 
                'electricity_price'    [0.2,0.2,0.2,0.2] 
                };
% Create DYNAMO parameters structure
params = DefaultOpts(params_in, default_params);


% Use poisson to get the probabilities of stochastic variables (outside temp, water
% usage demand, and non-controllable load demand)
params.stochastic_prob = cell(params.time_slots,1);
s_range_mat = cell2mat(params.stochastic_range);
 for t = 1:params.time_slots
    for p = 1:1
        s_min = s_range_mat(p,1);
        s_max = s_range_mat(p,2);
        s_unit = params.unit_space(p);
        params.stochastic_prob{t,p} = normpdf(s_min:s_unit:s_max,params.stochastic_medi(t,p),params.stochastic_sigma(p)*(params.stochastic_medi(t,p)+eps))';%Y = normpdf(X,MU,SIGMA),default(~,0,1)
        coef = 1/sum(params.stochastic_prob{t,p}); 
        params.stochastic_prob{t,p}=  params.stochastic_prob{t,p} *coef;% Sum of all probablities should be one.
    end
 end
end

