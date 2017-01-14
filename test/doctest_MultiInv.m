% Quick and dirty doctest style tests for MultiInv functions
%
% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   4  2016-10-27 15:01  BryanP      Added random sample related tests 
%   3  2016-10-27 13:16  BryanP      Added ops costs tests and corrected some others 
%   2  2016-10-27 12:46  BryanP      Overhauled to include both a runable script and the doctest comments 
%   1  2016-10           BryanP      Initial Code version or two 
%
% ============
%  DOC TESTS:
% ============
% Have to be part of first comment block
%
%% ===== Fix rand for repeatable results & Setup the problem =====
% >> rng('default')
% >> MultiInv_demo
% >> t = 1;
% >> test_pre_state = multi_inv_problem.state_set{t}.sample()
% 
% test_pre_state =
% 
%      8     5     2
% 
% 
% 
% % ===== Test Decision related functions ======
% >> dec_set_test = multi_inv_problem.fDecisionSet(multi_inv_problem.params, t, test_pre_state)
% 
% dec_set_test = 
% 
%   setBasic with properties:
% 
%     DStepsForContinuous: 10
%                  MinMax: [2x3 double]
%                StepSize: [1 1 1]
%                    name: ''
%            pt_dim_names: {}
%              SampleType: 'rand'
%                   N_dim: 3
%            DiscreteMask: [1 1 1]
% 
% >> dec_list = dec_set_test.as_array;
% >> size(dec_list)
% 
% ans =
% 
%    114     3
% 
% >> dec_subset = dec_list(1:10,:)
% 
% dec_subset =
% 
%      0     0     0
%      0     0     1
%      0     0     2
%      0     0     3
%      0     0     4
%      0     0     5
%      0     0     6
%      0     0     7
%      0     0     8
%      0     0     9
% 
% >> multi_inv_problem.fDecisionCost(multi_inv_problem.params, t, test_pre_state, dec_subset)
% 
% ans =
% 
%      0
%     -6
%     -8
%    -10
%    -12
%    -14
%    -16
%    -18
%    -20
%    -22
% 
% >> post_state_list = multi_inv_problem.fDecisionApply(multi_inv_problem.params, t, test_pre_state, dec_subset)
% 
% post_state_list =
% 
%      8     5     2
%      8     5     3
%      8     5     4
%      8     5     5
%      8     5     6
%      8     5     7
%      8     5     8
%      8     5     9
%      8     5    10
%      8     5    11
% 
% 
% 
% % ===== Test Uncertainty related functions ======
% >> sample_list = RandSetSample(multi_inv_problem.random_items, t, 10)
% 
% sample_list =
% 
%      2     2     6
%      1     3     4
%      0     1     4
%      0     2     4
%      1     4     2
%      3     3     4
%      3     5     1
%      0     2     4
%      3     0     0
%      3     3     2
% 
% >> rand_contrib = multi_inv_problem.fRandomCost(multi_inv_problem.params, t, post_state_list, sample_list)
% 
% rand_contrib =
% 
%     48
%     56
%     40
%     48
%     56
%     80
%     72
%     48
%     24
%     64
% 
% >> next_pre_list = multi_inv_problem.fRandomApply(multi_inv_problem.params, t, post_state_list, sample_list)
% 
% next_pre_list =
% 
%      6     3     0
%      7     2     0
%      8     4     0
%      8     3     1
%      7     1     4
%      5     2     3
%      5     0     7
%      8     3     5
%      5     5    10
%      5     2     9
%
%
% % ===== Test Ops cost Function in each of the 3 possible configurations =====
% % First check within multi_inv_problem (only 1 of 3 defined)
% >> isnan(multi_inv_problem.fOpsBeforeDecision)
% 
% ans =
% 
%      1
% 
% >> func2str(multi_inv_problem.fOpsAfterDecision)
% 
% ans =
% 
% MultiInvOps
% 
% >> isnan(multi_inv_problem.fOpsAfterRandom)
% 
% ans =
% 
%      1
% 
% 
% 
% % Now call underlying function in each instance to check behavior
% >> MultiInvOps(multi_inv_problem.params, t, test_pre_state)
% Warning: Unused Operations call before decision. Only after decision operations used in Multi-Inventory 
% *** 
% 
% ans =
% 
%      0
% 
% >> multi_inv_problem.fOpsAfterDecision(multi_inv_problem.params, t, post_state_list, dec_subset)
% 
% ans =
% 
%    -15
%    -16
%    -17
%    -18
%    -19
%    -20
%    -21
%    -22
%    -23
%    -24
% 
% >> MultiInvOps(multi_inv_problem.params, t, post_state_list, dec_subset, 'junk')
% Warning: Unused Operations after random. Only after decision operations used in Multi-Inventory 
% ***
% 
% ans =
% 
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
% 
% 
% 
% % And check full_state behavior
% >> [~, full_list] = MultiInvOps(multi_inv_problem.params, t, test_pre_state, [], [], 'old_full')
% Warning: Unused Operations call before decision. Only after decision operations used in Multi-Inventory 
% *** 
% 
% full_list =
% 
%      8     5     2
% 
% >> [~, full_list] = multi_inv_problem.fOpsAfterDecision(multi_inv_problem.params, t, post_state_list, dec_subset, [], 'old_full')
% 
% full_list =
% 
%      8     5     2
%      8     5     3
%      8     5     4
%      8     5     5
%      8     5     6
%      8     5     7
%      8     5     8
%      8     5     9
%      8     5    10
%      8     5    11
% 
% >> [~, full_list] = MultiInvOps(multi_inv_problem.params, t, post_state_list, dec_subset, 'junk', 'old_full')
% Warning: Unused Operations after random. Only after decision operations used in Multi-Inventory 
% *** 
% 
% full_list =
% 
%      8     5     2
%      8     5     3
%      8     5     4
%      8     5     5
%      8     5     6
%      8     5     7
%      8     5     8
%      8     5     9
%      8     5    10
%      8     5    11
% 
% 
% 
% % ===== Test Terminal Value Functions =====
% >> final_states = multi_inv_problem.state_set{multiinv_n_periods}.sample(10)
% 
% final_states =
% 
%      0     4     3
%      0     1     2
%      8     2     9
%      6     2    19
%      3     3     6
%      9     4    11
%      0     4     4
%      4     1    15
%      3     4     5
%      7     3    10
% 
% >> cache_term_val = multi_inv_problem.params.term_unit_val;
% >> multi_inv_problem.params.term_unit_val = 1:multi_inv_problem.params.n_products;
% >> MultiInvTermValue(multi_inv_problem.params, multi_inv_problem.n_periods, final_states)
% 
% ans =
% 
%     17
%      8
%     39
%     67
%     27
%     50
%     20
%     51
%     26
%     43
% 
% >> multi_inv_problem.params.term_unit_val = cache_term_val;
%







%% ============
% SCRIPT
% ============

%=== Script to run all commands once. Then can copy paste results into the
%following comments to run with doctest

% ===== Fix rand for repeatable results & Setup the problem =====
clear all
rng('default')
MultiInv_demo   %Runs MultiInvParamsSetup
echo on %Display commands as they are called for easy copy/paste. After demo to avoid excess output


t = 1;
test_pre_state = multi_inv_problem.state_set{t}.sample()


% ===== Test Decision related functions ======
dec_set_test = multi_inv_problem.fDecisionSet(multi_inv_problem.params, t, test_pre_state)
dec_list = dec_set_test.as_array;
size(dec_list)
dec_subset = dec_list(1:10,:)
multi_inv_problem.fDecisionCost(multi_inv_problem.params, t, test_pre_state, dec_subset)
post_state_list = multi_inv_problem.fDecisionApply(multi_inv_problem.params, t, test_pre_state, dec_subset)


% ===== Test Uncertainty related functions ======
sample_list = RandSetSample(multi_inv_problem.random_items, t, 10)
rand_contrib = multi_inv_problem.fRandomCost(multi_inv_problem.params, t, post_state_list, sample_list)
next_pre_list = multi_inv_problem.fRandomApply(multi_inv_problem.params, t, post_state_list, sample_list)


% ===== Test Ops cost Function in each of the 3 possible configurations =====
% First check within multi_inv_problem (only 1 of 3 defined)
isnan(multi_inv_problem.fOpsBeforeDecision)
func2str(multi_inv_problem.fOpsAfterDecision)
isnan(multi_inv_problem.fOpsAfterRandom)


% Now call underlying function in each instance to check behavior
MultiInvOps(multi_inv_problem.params, t, test_pre_state)
multi_inv_problem.fOpsAfterDecision(multi_inv_problem.params, t, post_state_list, dec_subset)
MultiInvOps(multi_inv_problem.params, t, post_state_list, dec_subset, 'junk')


% And check full_state behavior
[~, full_list] = MultiInvOps(multi_inv_problem.params, t, test_pre_state, [], [], 'old_full')
[~, full_list] = multi_inv_problem.fOpsAfterDecision(multi_inv_problem.params, t, post_state_list, dec_subset, [], 'old_full')
[~, full_list] = MultiInvOps(multi_inv_problem.params, t, post_state_list, dec_subset, 'junk', 'old_full')


% ===== Test Terminal Value Functions =====
final_states = multi_inv_problem.state_set{multiinv_n_periods}.sample(10)
cache_term_val = multi_inv_problem.params.term_unit_val;
multi_inv_problem.params.term_unit_val = 1:multi_inv_problem.params.n_products;
MultiInvTermValue(multi_inv_problem.params, multi_inv_problem.n_periods, final_states)
multi_inv_problem.params.term_unit_val = cache_term_val;


echo off  %Back to normal (to avoid printing all the doctest comments