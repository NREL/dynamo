function idx = MultiInvLookup(s_list, params)
% MULTIINVLOOKUP index (reverse) lookup from a given multi-d inventory state
%
% Usage: idx = MultiInvLookup(s_list, params)
%
% Returns a vector of state indexes corresponding the the multi-product
% inventory states in s_list. These numbers match those from
% MultiInvState().
%
% NOTE: May no longer be needed
%
% see also:
%   MultiInvState
%
% originally by Bryan Palmintier 2010

% IMPLEMENTATION NOTE:
% Note I tried multiple different approaches to this problem. All of these
% approaches achieve the same results, but take different amounts of time,
% as noted below. 

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2010-05-15 20:47  BryanP      original concept code
%   2  2010-05-16 01:45  BryanP      completed versions 1:3 and tested
%   3  2010-05-16 21:30  BryanP      Added #4: n-D matrix (fastest yet)


%% -- 1) reverse statelookup using either:
% rows = size(s_list,1);
% idx = zeros(rows,1);
% for s = 1:rows
%     % PERFORMANCE: Generating the key using sprintf takes up a large
%     % percentage of the time and is the root of the slow performance
%     key= sprintf('i%d',s_list(s,:));
%     % --     1a) a MATLAB structure
%     % PERFORMANCE During testing with random and orderd s_lists of all
%     % sizes, this method was consistantly the fastest
%     
%     idx(s) = params.state_lookup_struct.(key);
%     %% --     b) a java.util.Hashtable
%     % PERFORMANCE: Surprisingly this method takes more than 2.5x as long as
%     % the matlab structure lookup
%     %
%     % NOTE: for this sub-method to work, you must uncomment the java
%     % hashtable creation code in MultiInvInit
%     
%     %idx(s) = params.state_hash.get(key);
% end
%% -- 2) repeated use of find()
% PERFORMANCE: Surprisingly this method is fairly quick, but still takes
% about 2x as long as the MATLAB structure approach
%
% NOTE: the advantage of this techique is that it does not require storing
% a complete reverse lookup table.

% [rows,cols] = size(s_list);
% idx = zeros(rows,1);
% for r = 1:rows
%     match = ones(params.n_states,1);
%     for c = 1:cols;
%         match=and(match,(s_list(r,c) == params.states(:,c)));
%     end
%     idx(r) = find(match,1,'first');
% end

%% -- 3) using MATLABs intersect()
% IMPORTANT: this method only works if s_list contains no duplicates and if
% s_list is sorted in the same way as params.states
%
% PERFORMANCE: IF the above condition is met, this technique is the fastest
% method, clocking in about 2-2.5x faster than the MATLAB structure lookup.
%[junk, idx] = intersect(params.states, s_list,  'rows');

%% -- 4) Duh... just use an n-dimensional array
% PERFORMANCE: This is the fastest by far. It is about 3.5x faster than the
% MATLAB structure lookup. Since the inverse lookup is performed so often,
% this results in about a 20% speed-up of the over-all multi-dimensional
% DP!
%
% The triky part is how to access an array of subscripts with out looping.
% Thankfully, MATLAB provides sub2ind(), which with a little work lets us
% do just that.

% Build up a cell array containing the columns of the states to lookup
%
% Note: add one b/c the 0 inventory level is an invalid MATLAB subscript
subscript_list = cell(1,params.num_products);
for p = 1:params.num_products
    subscript_list{p} = s_list(:,p) + 1;
end

% And now extract the equivalent single index within the sparce state
% lookup table for all of our states.
%
% Note: using {:} effectively converts each element of subscript_list
% to a separate parameter, thereby providing one column for each dimension
% as is required by sub2ind.
idx_list = sub2ind(size(params.state_lookup), subscript_list{:});

% finally extract the value stored at the corresponding single index, which
% was setup in MultiInvInit() to be the corresponding state number
idx = params.state_lookup(idx_list);
