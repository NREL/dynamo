function [next_pre_state_list] = ...
    SimpleStoragePvRandomApply(params, t, post_state_list, uncertain_pv_list) %#ok<INUSL>
% SIMPLESTORAGEPVRANDOMAPPLY Applies the provided uncertainty(ies) to the current (post-decision) state 
% 
% Produces a (set of) new pre decision state(s) to be used during the next time period. 

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-07-18 07:09  BryanP      Initial Code adapted from MultiInvRandomApply v1 

%Start with prost decision list, since much will be unchanged
next_pre_state_list = post_state_list;
if size(post_state_list,1) == 1
    next_pre_state_list = repmat(next_pre_state_list, size(uncertain_pv_list, 1),1);
end
% Swap in new pv levels
%   Note: the pv uncertainty is not additions, but rather the new
%   penetration level
next_pre_state_list(:,1) = uncertain_pv_list;
% And prevent any negative quantities
next_pre_state_list = max(0, next_pre_state_list);
