% Quick and dirty doctest style tests for MultiInv functions
%
% Note: test for runnind dpBI with thes included in the MultiInv_demo file
% as its own doctest
%
% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   7  2017-06-01 22:15  BryanP      Updated for SetCombinWithLimits based functions 
%   6  2017-04-03 10:25  BryanP      Corrected tests for other recent changes 
%   5  2017-04-03 10:15  BryanP      Updated to call MultiInv_demo as a function  
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
% >> multi_inv_problem = MultiInv_demo('med');   %Runs MultiInvParamsSetup
% 
% 
% >> t = 1;
% >> test_pre_state_list = multi_inv_problem.state_set{t}.sample(3)
% 
% test_pre_state_list =
% 
%      4     0    12
%      3     1     6
%      6     0     0
% 
% >> test_pre_state = multi_inv_problem.state_set{t}.sample()
% 
% test_pre_state =
% 
%      2     2     2
% 
% 
% 
% %% ===== Test Decision related functions ======
% >> dec_set_test = multi_inv_problem.fDecisionSet(multi_inv_problem.params, t, test_pre_state)
% 
% dec_set_test = 
% 
%   setCombinWithLimits with properties:
% 
%     UseValueSpace: 0
%               Opt: [1×1 struct]
%          N_Combin: 41
%              name: ''
%      pt_dim_names: {}
%        SampleType: 'from_list'
%             N_dim: 3
%      DiscreteMask: []
% 
% 
% >> dec_list = dec_set_test.as_array;
% >> size(dec_list)
% 
% ans =
% 
%     41     3
% 
% >> dec_subset = dec_list(1:10,:)
% 
% dec_subset =
% 
%      0     0     0
%      0     1     0
%      0     2     0
%      1     0     0
%      1     1     0
%      1     2     0
%      2     0     0
%      2     1     0
%      3     0     0
%      4     0     0
% 
% >> multi_inv_problem.fDecisionCost(multi_inv_problem.params, t, test_pre_state, dec_subset)
% 
% ans =
% 
%      0
%     -6
%     -8
%     -6
%     -8
%    -10
%     -8
%    -10
%    -10
%    -12
% 
% >> post_state_list = multi_inv_problem.fDecisionApply(multi_inv_problem.params, t, test_pre_state, dec_subset)
% 
% post_state_list =
% 
%      2     2     2
%      2     3     2
%      2     4     2
%      3     2     2
%      3     3     2
%      3     4     2
%      4     2     2
%      4     3     2
%      5     2     2
%      6     2     2
% 
% 
% 
% %% ===== Test Uncertainty related functions ======
% >> sample_list = RandSetSample(multi_inv_problem.random_items, 10, t)
% 
% sample_list =
% 
%      1     4     2
%      3     3     4
%      3     5     1
%      0     2     4
%      3     0     0
%      3     3     2
%      1     4     0
%      2     3     1
%      0     3     5
%      1     3     4
% 
% >> rand_contrib = multi_inv_problem.fRandomCost(multi_inv_problem.params, t, post_state_list, sample_list)
% 
% rand_contrib =
% 
%     40
%     56
%     56
%     32
%     24
%     64
%     24
%     48
%     32
%     40
% 
% >> next_pre_list = multi_inv_problem.fRandomApply(multi_inv_problem.params, t, post_state_list, sample_list)
% 
% next_pre_list =
% 
%      1     0     0
%      0     0     0
%      0     0     1
%      3     0     0
%      0     3     2
%      0     1     0
%      3     0     2
%      2     0     1
%      5     0     0
%      5     0     0
% 
% 
% 
% %% ===== Test Ops cost Function in each of the 3 possible configurations =====
% % First check within multi_inv_problem (only 1 of 3 defined)
% >> isempty(multi_inv_problem.fOpsBeforeDecision)
% 
% ans =
% 
%   logical
% 
%    1
% 
% >> func2str(multi_inv_problem.fOpsAfterDecision)
% 
% ans =
% 
% 'MultiInvOps'
% 
% >> isempty(multi_inv_problem.fOpsAfterRandom)
% 
% ans =
% 
%   logical
% 
%    1
% 
% 
% 
% %% Now call underlying function in each instance to check behavior
% >> MultiInvOps(multi_inv_problem.params, t, test_pre_state)
% Warning: Unused Operations call before decision. Only after decision operations used in Multi-Inventory 
% > In MultiInvOps (line 26)
%   **** 
% 
% >> multi_inv_problem.fOpsAfterDecision(multi_inv_problem.params, t, post_state_list, dec_subset)
% 
% ans =
% 
%     -6
%     -7
%     -8
%     -7
%     -8
%     -9
%     -8
%     -9
%     -9
%    -10
% 
% >> MultiInvOps(multi_inv_problem.params, t, post_state_list, dec_subset, 'junk')
% Warning: Unused Operations after random. Only after decision operations used in Multi-Inventory 
% > In MultiInvOps (line 33)
% ****
%
% 
% %% And check full_state behavior
% >> [~, full_list] = MultiInvOps(multi_inv_problem.params, t, test_pre_state, [], [], 'old_full')
% Warning: Unused Operations call before decision. Only after decision operations used in Multi-Inventory 
% > In MultiInvOps (line 26)
%   ****
% 
% >> [~, full_list] = multi_inv_problem.fOpsAfterDecision(multi_inv_problem.params, t, post_state_list, dec_subset, [], 'old_full')
% 
% full_list =
% 
%      2     2     2
%      2     3     2
%      2     4     2
%      3     2     2
%      3     3     2
%      3     4     2
%      4     2     2
%      4     3     2
%      5     2     2
%      6     2     2
% 
% >> [~, full_list] = MultiInvOps(multi_inv_problem.params, t, post_state_list, dec_subset, 'junk', 'old_full')
% Warning: Unused Operations after random. Only after decision operations used in Multi-Inventory 
% > In MultiInvOps (line 33)
%   **** 
% 
% 
% %% ===== Test Terminal Value Functions =====
% >> final_states = multi_inv_problem.state_set{multi_inv_problem.n_periods}.sample(10)
% 
% final_states =
% 
%      5     2     2
%      2     0    14
%      1     5     0
%      0     3     4
%      2     3     3
%      6     0     8
%      2     1     9
%      3     3     1
%      4     1     4
%      0     5     4
% 
% >> cache_term_val = multi_inv_problem.params.term_unit_val;
% >> multi_inv_problem.params.term_unit_val = 1:multi_inv_problem.params.n_products;
% >> MultiInvTerminalValue(multi_inv_problem.params, multi_inv_problem.n_periods, final_states)
% 
% ans =
% 
%     15
%     44
%     11
%     18
%     17
%     30
%     31
%     12
%     18
%     22
% 
% >> multi_inv_problem.params.term_unit_val = cache_term_val;







%% ============
% SCRIPT
% ============
%=== Script to run all commands once. Then can copy paste results into the
%following comments to run with doctest

%% ===== Fix rand for repeatable results & Setup the problem =====
clear all
echo on %Display commands as they are called for easy copy/paste. After demo to avoid excess output
rng('default')
multi_inv_problem = MultiInv_demo('med');   %Runs MultiInvParamsSetup


t = 1;
test_pre_state_list = multi_inv_problem.state_set{t}.sample(3)
test_pre_state = multi_inv_problem.state_set{t}.sample()


%% ===== Test Decision related functions ======
dec_set_test = multi_inv_problem.fDecisionSet(multi_inv_problem.params, t, test_pre_state)
dec_list = dec_set_test.as_array;
size(dec_list)
dec_subset = dec_list(1:10,:)
multi_inv_problem.fDecisionCost(multi_inv_problem.params, t, test_pre_state, dec_subset)
post_state_list = multi_inv_problem.fDecisionApply(multi_inv_problem.params, t, test_pre_state, dec_subset)


%% ===== Test Uncertainty related functions ======
sample_list = RandSetSample(multi_inv_problem.random_items, 10, t)
rand_contrib = multi_inv_problem.fRandomCost(multi_inv_problem.params, t, post_state_list, sample_list)
next_pre_list = multi_inv_problem.fRandomApply(multi_inv_problem.params, t, post_state_list, sample_list)


%% ===== Test Ops cost Function in each of the 3 possible configurations =====
% First check within multi_inv_problem (only 1 of 3 defined)
isempty(multi_inv_problem.fOpsBeforeDecision)
func2str(multi_inv_problem.fOpsAfterDecision)
isempty(multi_inv_problem.fOpsAfterRandom)


%% Now call underlying function in each instance to check behavior
MultiInvOps(multi_inv_problem.params, t, test_pre_state)
multi_inv_problem.fOpsAfterDecision(multi_inv_problem.params, t, post_state_list, dec_subset)
MultiInvOps(multi_inv_problem.params, t, post_state_list, dec_subset, 'junk')


%% And check full_state behavior
[~, full_list] = MultiInvOps(multi_inv_problem.params, t, test_pre_state, [], [], 'old_full')
[~, full_list] = multi_inv_problem.fOpsAfterDecision(multi_inv_problem.params, t, post_state_list, dec_subset, [], 'old_full')
[~, full_list] = MultiInvOps(multi_inv_problem.params, t, post_state_list, dec_subset, 'junk', 'old_full')


%% ===== Test Terminal Value Functions =====
final_states = multi_inv_problem.state_set{multi_inv_problem.n_periods}.sample(10)
cache_term_val = multi_inv_problem.params.term_unit_val;
multi_inv_problem.params.term_unit_val = 1:multi_inv_problem.params.n_products;
MultiInvTerminalValue(multi_inv_problem.params, multi_inv_problem.n_periods, final_states)
multi_inv_problem.params.term_unit_val = cache_term_val;


echo off  %Back to normal (to avoid printing all the doctest comments