function [ HEMS_problem,results,time ] = HEMS_Main( varargin )
%function [ problem,results,time] = HEMS_Main( bool,time_slot,initial_state,medi,price,ad_time,state_sample,decision_sample,uncertain_sample )
%Adapted from MultiInv by BryanP
%Changed by Li Wang
HEMS_discount = 0;    % Factor applied for future discounted values. 
HEMS_timeslots = 24;  % Scheduling horizon.
initial_state = [80,110,52]; % Default initial state
HEMS_medi =0; % Default mean for normal distruibution
HEMS_price = 0; % Default electricity price
HEMS_adtime = [6,20]; % Default arrive departure time of EV
if (~any(strcmpi('default', varargin)))
    if(~isnan(varargin{2}))
        HEMS_timeslots = varargin{2};
    end
    
    if(sum(isnan(varargin{3}))==0)
        initial_state = varargin{3};
    end
    
    if(sum(sum(isnan(varargin{4})))==0)
        HEMS_medi = varargin{4};
    end
    
    if(sum(isnan(varargin{5}))==0)
        HEMS_price = varargin{5};
    end
    
    if(sum(isnan(varargin{6}))==0)
        HEMS_adtime = varargin{6};
    end
end


% Home Energy Management System(HEMS) parameters.
HEMS_params = {
                'time_slots'              HEMS_timeslots
                'total_space'             realmax('single') %Big number to show there is no constraint on the total space.
                'unit_space'              [0.1] % Shows the amount of space ( degree temp., degree temp., %state of charge ) each unit occupies. 
                %HEMS_DecisionApply RandomApply round hardcode.
                'appliance_range'         {[70 85]} % The range using which state space for HVAC (air temp),WH (hot water temp),EV(SOC)
                                                                       %is generated.
                'stochastic_range'        {[80,110]} % The range of stochastic variables (outside temp, hot water usage, non-controllable load).
                                                                   %The non-controllable load unit is converted to percent of EV SOC.
                'stochastic_medi'          HEMS_medi% This is the [#of timeslots , #of appliances] matrix that shows the time series for the
                                                                       %most frequently realized stochastic variable. 
                'stochastic_sigma'        [0.1]
                 'unit_cost'              [-5.3] % The electricity usage that shows the of amount of kW needed to change the air/water temprature
                                                                %or SOC by one unit.
                 'desired_value'          [75] % Desired values for the appliances.
                 %'arrive_departure_time'  HEMS_adtime % Arrival and departure time of EV.
                 %'electricity_price'     [0.1,0.1,0.1,0.1,1,1,1]
                 'electricity_price'      HEMS_price % The [1, #of timeslots] electricity time series price for every time slot. 
                };
           


 
state_sample = 500;
decision_sample =30;
uncertain_sample = 30;  

    if(nargin>6&&~isnan(varargin{7}))
        state_sample = varargin{7};
    end

    if(nargin>7&&~isnan(varargin{8}))
        decision_sample = varargin{8};
    end

    if(nargin>8&&~isnan(varargin{9}))
       uncertain_sample = varargin{9};
    end

 % These parameters are used by the adpSBI function.
sbi_opt = {     'sbi_state_samples_per_time'            state_sample    % Number of state samples per time period
                'sbi_decisions_per_sample'              decision_sample    % Number of decision samples per state
                'sbi_uncertain_samples_per_post'        uncertain_sample  % Number of random/uncertainty samples per time, used for all decisions
                'vfun_approx'                           'LocalRegr'%LocalRegr,LocalAvg
                };
            


HEMS_params = HEMS_setup_params(HEMS_params); % Set up HEMS parameters.
HEMS_problem = HEMS_setup_problem(HEMS_params,'initial_state',initial_state,'discount_rate',HEMS_discount,'n_periods',HEMS_timeslots); % Set up HEMS problem

% To run dynamic programming using backward induction (dpBI)
if (nargin>0 && islogical(varargin{1})&& varargin{1} == true)
    tic;
    results = dpBI(HEMS_problem);
    time = toc;

% To run adaptive dynamic programming using sampled backward
% induction(adpSBI)
else
tic;
%results = dpBI(HEMS_problem);
results = adpSBI(HEMS_problem, sbi_opt);
time = toc;
end
end

