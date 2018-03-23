function [simple_storage_pv_problem, results] = SimpleStoragePv_demo(varargin)
% SIMPLESTORAGEPV_DEMO Demo script showing simple use of adp clases of algorithms
%
% simple_storage_pv_problem = SimpleStoragePv_demo()
%   set up full example simple_storage_pv problem structure
%
% SimpleStoragePv_demo(scenario_str)
%   specify which scenario to setup. Valid options are:
%     'isgt'         small, problem that matches Palmintier, et al 2016 using on-line calls to OpenDSS 
%     'isgt_lookup'  same as isgt, but uses pre-computed lookup table of OpenDSS output results (Default) 
%     'medium'       A larger problem 
%
% SimpleStoragePv_demo(__, 'dp')
% SimpleStoragePv_demo(__, 'dpBI')
%   Run traditional backward induction on example SimpleStoragePv problem
%
% SimpleStoragePv_demo(__, 'sbi')
% SimpleStoragePv_demo(__, 'adpSBI')
%   Run sampled backward induction ADP algorithm on example SimpleStoragePv problem
%
% [simple_storage_pv_problem, results] = SimpleStoragePv_demo(__)
%   Return the result of the corresponding simulation as a variable
%
% Example and doctests:
%
% >> format shortG
% >> simple_storage_pv_problem = SimpleStoragePv_demo()
% 
% simple_storage_pv_problem = 
%  struct with fields:
% 
%                     params: [1×1 struct]
%              discount_rate: 0.15065
%                  n_periods: 2
%                  state_set: {[1×1 setList]  [1×1 setList]  [1×1 setList]}
%             fTerminalValue: @SimpleStoragePvTerminalValue
%               fDecisionSet: @SimpleStoragePvDecisionSet
%              fDecisionCost: @SimpleStoragePvDecisionValue
%             fDecisionApply: @SimpleStoragePvDecisionApply
%          decision_vfun_map: []
%               random_items: {[1×1 rpTransMatrix]}
%                fRandomCost: @SimpleStoragePvRandomValue
%               fRandomApply: @SimpleStoragePvRandomApply
%           random_state_map: {[1]}
%         fOpsBeforeDecision: @SimpleStoragePvOps
%          fOpsAfterDecision: @SimpleStoragePvOps
%            fOpsAfterRandom: []
%     fCompareDecisionPolicy: []
%              fCompareValue: []
%             fMapState2Vfun: []
%           fOptimalDecision: []
%              fRandomSample: []
%               fRandomJoint: []
%
% >> [~, isgt2016_results] = SimpleStoragePv_demo('isgt_lookup', 'dp');
% Backward Induction DP
%     T=3 (terminal period): Done
%     T=2:Done: 6 states
%     T=1:Done: 1 states
% Elapsed time is *** seconds.
%
% >> isgt2016_results
% 
% isgt2016_results = 
% 
%   struct with fields:
% 
%        dpbi_policy: {[1000 1000]  [6×2 double]}
%            dp_opts: [1×1 struct]
%        dpbi_values: {[-1.0804e+06]  [6×1 double]  [9×1 double]}
%     pre_state_list: {[0.76 0 0]  [6×3 double]  [9×3 double]}
% 
% >> [~ , isgt2016_as_adp] = SimpleStoragePv_demo('isgt_lookup', 'sbi');
% Sampled Backward Induction ([0.1] ksamples/period)
%     Creating empty post-decision value functions (LocalAvg)
%     T=3 (terminal period): Done
%     T=2:S........................................100
%     T=1:S........................................100
% Elapsed time is ***
% 
% >> isgt2016_as_adp
%  isgt2016_as_adp = 
% 
%   struct with fields:
% 
%                log: [1×1 struct]
%          post_vfun: [1×3 faLocalAvg]
%            adp_opt: [1×1 struct]
%     first_decision: [1000 1000]
%          objective: ***
%
% >> relative_abs_error=abs(isgt2016_as_adp.objective +  1080352)/1080352
% relative_abs_error =
%    ***
% 
% relative_abs_error < 0.01
% 
% ans =
%   logical
%    1
%
%
% ISGT Reference:
%  B. Palmintier, D. Krishnamurthy, and H. Wu, ?Design Flexibility for
%  Uncertain Distributed Generation from Photovoltaics,? in Innovative
%  Smart Grid Technologies Conference 2016, Minneapolis, MN, 2016.


% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   4  2018-03-23 02:36  BryanP      Default to using lookup table, additional comments, doctest isgt_lookup results 
%   3  2017-07-17 17:34  BryanP      Added ISGT lookup data table 
%   2  2017-07-16 21:24  BryanP      Completed initial parameter, etc. definition 
%   1  2017-07-11 19:07  BryanP      Adapted from MultiInv_demo v1

%% Setup problem specific values
if any(strcmpi('medium', varargin)) || any(strcmpi('med', varargin))
    % Create "medium" test case if requested.
    error('SimpStorPv:not_implemented','Medium Storage Sizing for PV problem not yet implemented')

elseif any(strcmpi('isgt_lookup', varargin))
    % ISGT-2016 PROBLEM using PRE-COMPUTED LOOKUP Operations

    SimpleStoragePv_discount_per_yr = 0.04;
	%Use pre-computed operations from ISG problem
    stored_ops = [                      %Net Energy in kWh 
      % PV_%      Bat_kW    Bat       Bckfd-kWh Hour 1    Hour 2    Hour 3    Hour 4    Hour 5    Hour 6    Hour 7    Hour 8    Hour 9    Hour 10   Hour 11   Hour 12   Hour 13   Hour 14   Hour 15   Hour 16   Hour 17   Hour 18   Hour 19   Hour 20   Hour 21   Hour 22   Hour 23   Hour 24
        0.76      0         0         0.000     4062.336  3944.626  3849.338  3916.779  4202.711  4884.138  5745.423  6405.025  5736.520  4974.035  3816.175  3164.878  3437.039  4003.031  4475.898  4568.945  5681.599  6838.350  6519.956  6160.491  5758.631  5211.111  4638.089  4259.734
        0.76      200       1000      0.000     4435.690  4344.628  4303.627  4370.990  4495.452  4846.803  5404.198  5827.382  5344.865  4839.051  4249.168  3905.790  4004.612  4370.169  4617.488  4669.991  5296.307  6080.231  5757.624  5484.204  5309.784  5094.268  4664.098  4520.294
        0.76      1000      1000      0.000     4354.747  4263.894  4202.200  4269.605  4495.507  4884.461  5513.924  5723.545  5280.930  4865.664  4186.946  4018.411  4117.104  4452.567  4547.118  4587.068  5390.158  5835.740  5753.414  5684.302  5505.499  5166.397  4664.090  4468.807
        0.76      1000      5000      0.000     4852.403  4816.725  4809.795  4834.752  4862.217  4857.833  4825.961  4861.185  4850.410  4821.102  4826.359  4780.855  4771.889  4827.104  4845.576  4841.364  4854.302  4907.485  4871.705  4847.897  4815.069  4827.689  4858.062  4867.451
        0.86      200       1000      0.000     4315.395  4242.897  4147.240  4214.652  4430.112  4884.458  5329.875  5747.067  5076.425  4444.144  3800.379  3197.869  3274.425  3808.995  4056.794  4106.427  5002.116  6027.416  5706.659  5507.356  5294.429  5068.404  4683.845  4410.424
        0.86      1000      1000      0.000     4273.539  4252.340  4156.726  4224.136  4388.412  4884.456  5454.833  5626.551  5040.051  4407.892  3683.064  3427.586  3435.053  3960.904  3986.225  4014.069  5077.578  5714.867  5658.051  5573.773  5436.096  5151.220  4637.653  4368.181
        0.86      1000      5000      0.000     4630.886  4641.873  4612.841  4630.776  4628.438  4615.820  4593.200  4640.722  4631.246  4622.144  4620.622  4506.613  4504.329  4610.462  4632.387  4623.179  4620.037  4656.059  4651.712  4608.728  4579.793  4605.372  4626.599  4638.713
        1.09      200       1000      31.204    4127.427  4022.594  3962.073  4009.034  4267.274  4849.411  5188.713  5584.346  4629.496  3734.705  2398.831  1465.521  1786.272  2269.154  2897.518  3024.914  4346.637  5867.974  5677.244  5441.890  5259.652  5073.945  4651.716  4324.400
        1.09      1000      1000      0.000     4062.378  3944.626  3849.338  3916.779  4202.711  4884.138  5361.452  5558.393  4591.097  3536.668  2294.965  2090.928  2028.805  2471.361  2672.190  2760.807  4397.837  5606.823  5629.365  5530.133  5409.243  5114.482  4637.652  4259.323
        1.09      1000      5000      0.000     4139.966  4137.398  4158.802  4138.802  4143.353  4103.769  4116.614  4163.983  4084.878  4089.658  4034.450  3911.576  3863.969  4026.847  4046.907  4057.173  4098.990  4144.315  4218.369  4164.277  4096.948  4091.785  4120.825  4130.919
        1.11      200       1000      94.214    4139.965  4022.312  3941.283  4008.745  4279.899  4849.412  5185.877  5551.894  4601.107  3657.739  2284.179  1339.906  1664.833  2144.871  2791.564  2934.994  4288.153  5849.728  5651.765  5454.939  5273.456  5073.946  4630.925  4338.644
        1.11      1000      1000      0.000     4062.378  3944.626  3849.338  3916.779  4202.711  4884.138  5358.617  5478.324  4561.813  3462.766  2188.369  1958.908  1900.375  2351.245  2563.982  2668.306  4341.030  5622.855  5607.331  5546.619  5426.484  5114.483  4637.652  4259.323
        1.11      1000      5000      0.000     4110.391  4115.984  4118.001  4115.236  4114.635  4087.340  4090.182  4088.097  3997.983  4012.420  3996.631  3888.880  3834.361  3991.345  4005.564  4034.264  4054.021  4115.460  4168.287  4142.955  4081.736  4075.172  4101.381  4113.439
        1.26      200       1000      1572.797  4062.378  3944.626  3849.338  3916.779  4202.711  4849.878  5180.748  5565.487  4266.104  2971.501  1323.775  145.425   505.184   1111.865  1950.533  2101.304  3827.024  5713.628  5608.408  5438.045  5318.844  5057.086  4603.376  4259.322
        1.26      1000      1000      1260.737  4062.378  3944.626  3849.338  3916.779  4202.711  4884.138  5249.317  5398.692  4189.792  2789.135  1089.954  850.426   778.315   1334.491  1565.027  1816.982  3871.172  5594.645  5535.730  5451.774  5375.956  5103.691  4637.651  4259.323
        1.26      1000      5000      0.000     3918.735  3921.025  3915.459  3915.772  3921.443  3903.022  3892.901  3991.981  3647.466  3352.133  3327.282  3092.906  3080.705  3305.209  3384.551  3390.502  3580.765  4044.356  4023.586  3958.869  3903.443  3903.078  3919.236  3916.313
        1.26      200       1000      1579.605  4062.378  3944.626  3849.338  3916.779  4202.711  4849.878  5180.669  5565.023  4262.731  2967.545  1308.975  124.413   489.188   1098.351  1935.716  2092.168  3821.478  5712.842  5608.325  5438.021  5318.844  5057.086  4603.376  4259.322
        1.26      1000      5000      0.000     3911.010  3916.643  3920.878  3923.387  3921.610  3894.593  3883.313  3990.374  3633.056  3344.175  3307.163  3079.486  3060.431  3295.487  3363.343  3374.384  3574.899  4042.892  4021.996  3957.309  3894.953  3901.524  3910.826  3908.475
        1.26      1000      5000      0.000     3911.010  3916.643  3920.878  3923.387  3921.610  3894.593  3883.313  3990.374  3633.056  3344.175  3307.163  3079.486  3060.431  3295.487  3363.343  3374.384  3574.899  4042.892  4021.996  3957.309  3894.953  3901.524  3910.826  3908.475
        1.51      200       1000      7585.338  4062.378  3944.626  3849.338  3916.779  4202.711  4843.019  5090.532  5450.921  3641.231  1785.378  -354.746  -1862.281 -1452.466 -737.624  293.469   661.387   2950.848  5463.465  5525.594  5406.094  5199.986  5034.556  4596.514  4259.322
        1.51      1000      5000      664.115   3724.101  3732.884  3704.394  3714.526  3721.690  3701.490  3678.028  3874.989  3194.724  2335.876  1973.678  1278.419  1485.880  1861.379  2036.179  2029.951  2676.843  3840.096  3816.905  3776.665  3724.502  3706.108  3714.453  3721.264
        1.51      1000      5000      664.115   3724.101  3732.884  3704.394  3714.526  3721.690  3701.490  3678.028  3874.989  3194.724  2335.876  1973.678  1278.419  1485.880  1861.379  2036.179  2029.951  2676.843  3840.096  3816.905  3776.665  3724.502  3706.108  3714.453  3721.264
        1.51      200       1000      7585.338  4062.378  3944.626  3849.338  3916.779  4202.711  4843.019  5090.532  5450.921  3641.231  1785.378  -354.746  -1862.281 -1452.466 -737.624  293.469   661.387   2950.848  5463.465  5525.594  5406.094  5199.986  5034.556  4596.514  4259.322
        1.51      1000      1000      7042.513  4062.378  3944.626  3849.338  3916.779  4202.711  4869.462  5166.767  5296.816  3558.995  1639.169  -766.281  -996.346  -1117.383 -556.497  -145.171  359.112   3052.274  5308.903  5439.579  5380.676  5295.515  5092.401  4637.651  4259.323
        1.51      1000      5000      664.115   3724.101  3732.884  3704.394  3714.526  3721.690  3701.490  3678.028  3874.989  3194.724  2335.876  1973.678  1278.419  1485.880  1861.379  2036.179  2029.951  2676.843  3840.096  3816.905  3776.665  3724.502  3706.108  3714.453  3721.264
        0.00      0         0         0.000     4062.754  3945.029  3849.317  3916.779  4202.711  4884.138  5870.452  6931.723  7836.904  8376.839  8828.875  9089.869  9214.764  9432.405  9316.700  8868.160  8383.358  7664.918  6742.236  6202.417  5762.774  5212.317  4639.250  4260.159
        ];
    SimpleStoragePv_params = { %REQUIRED: n_periods, storage_kVA, storage_kWh
                        'n_periods'         2                     % Number of decision periods (does not include terminal period
                        'storage_kVA'       {0, [ 200; 1000]}     % List of valid storage power ratings (kVA) as a function of time
                        'storage_kWh'       {0, [1000; 5000]}     % List of valid storage capacities (kWh)
                        'storage_hr_limits' [1, 5]                % Limit storage power/energy combinations based on hour rating

                        'run_opendss'       false       % Use pre-computed results
                        'stored_ops'        stored_ops  % Table of pre-computed values 
                        
                        'scrap_value_on_inv_upsize' [0 0]   % Map of storage components with scrap value if inverter increases in size
                    };

    sbi_opt = { 'sbi_state_samples_per_time'           100    % Number of state samples per time period
                'sbi_decisions_per_sample'              20     % Number of decision samples per state
                'sbi_uncertain_samples_per_post'        10     % Number of random/uncertainty samples per time, used for all decisions
                
                'vfun_approx'                           'LocalAvg'
                };
            
elseif any(strcmpi('isgt', varargin))
    % ISGT-2016 PROBLEM with on-line OpenDSS runs

    SimpleStoragePv_discount_per_yr = 0.04;
    SimpleStoragePv_params = { %REQUIRED: n_periods, storage_kVA, storage_kWh
                        'n_periods'         2                     % Number of decision periods (does not include terminal period
                        'storage_kVA'       {0, [ 200; 1000]}     % List of valid storage power ratings (kVA) as a function of time
                        'storage_kWh'       {0, [1000; 5000]}     % List of valid storage capacities (kWh)
                        'storage_hr_limits' [1, 5]                % Limit storage power/energy combinations based on hour rating
                        'run_opendss'       true        % Run OpenDSS for results
                    };
                 
    sbi_opt = { 'sbi_state_samples_per_time'           100    % Number of state samples per time period
                'sbi_decisions_per_sample'              20     % Number of decision samples per state
                'sbi_uncertain_samples_per_post'        10     % Number of random/uncertainty samples per time, used for all decisions
                
                'vfun_approx'                           'LocalAvg'
                };
else %No known problem definition
    %Default to using isgt_lookup
    if nargout < 2
        simple_storage_pv_problem = SimpleStoragePv_demo('isgt_lookup', varargin{:});
    else
        [simple_storage_pv_problem, results] = SimpleStoragePv_demo('isgt_lookup', varargin{:});
    end
    return
end

% Build up additional derived fields. 
SimpleStoragePv_params = SimpleStoragePvSetupParams(SimpleStoragePv_params);
yrs_per_period = SimpleStoragePv_params.ops_pre_decision_yrs + SimpleStoragePv_params.ops_post_decision_yrs;
simple_storage_pv_problem = SimpleStoragePvSetupProblem(SimpleStoragePv_params, ...
    'discount_rate', 1-(1-SimpleStoragePv_discount_per_yr)^yrs_per_period, 'n_periods', SimpleStoragePv_params.n_periods);
             
%% Now that setup is complete, let's run the specified examples
if nargin < 1
    return
end


%% Run DP (Backward Induction) algorithm
if any(strcmpi('dp', varargin)) || any(strcmpi('dpBI', varargin))
    simple_storage_pv_problem_dp = simple_storage_pv_problem;

    tic
    results = dpBI(simple_storage_pv_problem_dp);
    toc
    %TODO: display summary
    return
end

%% Run ADP Sample Backward Induction algorithm
if any(strcmpi('sbi', varargin)) || any(strcmpi('adpSBI', varargin))
    simple_storage_pv_problem_sbi = simple_storage_pv_problem;

    tic
    results = adpSBI(simple_storage_pv_problem_sbi, sbi_opt);
    toc
    %TODO: display summary
    return
end

%% Run ADP Double Pass (Temporal Difference, Lambda=1) algorithm
if any(strcmpi('td1', varargin)) || any(strcmpi('adpTD1', varargin))
    simple_storage_pv_problem_td1 = simple_storage_pv_problem;
    tic
    results = adpTD1(simple_storage_pv_problem_td1);
    toc
    %TODO: display summary
    return
end

