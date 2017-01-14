% Quick & dirty demonstration of RandPartition
%
% consider stepping through with debugger

%% Simple 4-D calls for visual inspection
p = RandPartition([2 3 5 7], 30*ones(10,1),'lumpy',false), disp(nnz(p))
p = RandPartition([2 3 5 7], 30*ones(10,1),'lumpy',true), disp(nnz(p))
p = RandPartition([2 3 5 7], 30*ones(10,1),'lumpy',true,'lump_factor',1), disp(nnz(p))
p = RandPartition([2 3 5 7], 30*ones(10,1),'lumpy',true,'lump_factor',2), disp(nnz(p))
p = RandPartition([2 3 5 7], 30*ones(10,1),'lumpy',true,'lump_factor',3), disp(nnz(p))
p = RandPartition([2 3 5 7], 30*ones(10,1),'lumpy',true,'lump_factor',4), disp(nnz(p))

%% 6-D calls, visual still OK
p = RandPartition([2 3 5 7 11 13], 100*ones(10,1),'lumpy',false), disp(nnz(p))
p = RandPartition([2 3 5 7 11 13], 100*ones(10,1),'lumpy',true), disp(nnz(p))
p = RandPartition([2 3 5 7 11 13], 100*ones(10,1),'lumpy',true,'lump_factor',1), disp(nnz(p))
p = RandPartition([2 3 5 7 11 13], 100*ones(10,1),'lumpy',true,'lump_factor',2), disp(nnz(p))
p = RandPartition([2 3 5 7 11 13], 100*ones(10,1),'lumpy',true,'lump_factor',3), disp(nnz(p))
p = RandPartition([2 3 5 7 11 13], 100*ones(10,1),'lumpy',true,'lump_factor',4), disp(nnz(p))
p = RandPartition([2 3 5 7 11 13], 100*ones(10,1),'lumpy',true,'lump_factor',5), disp(nnz(p))
p = RandPartition([2 3 5 7 11 13], 100*ones(10,1),'lumpy',true,'lump_factor',6), disp(nnz(p))

%% 6-D Larger sample sizes (raw results hidden) for lumpiness comparisions
n=1000; p = RandPartition([2 3 5 7 11 13], 100*ones(n,1),'lumpy',true,'lump_factor',[]); disp(nnz(p)/n)

n=1000; p = RandPartition([2 3 5 7 11 13], 100*ones(n,1),'lumpy',true,'lump_factor',1); disp(nnz(p)/n)
n=1000; p = RandPartition([2 3 5 7 11 13], 100*ones(n,1),'lumpy',true,'lump_factor',2); disp(nnz(p)/n)
n=1000; p = RandPartition([2 3 5 7 11 13], 100*ones(n,1),'lumpy',true,'lump_factor',6); disp(nnz(p)/n)

n=1000; p = RandPartition([2 3 5 7 11 13], 100*ones(n,1),'lumpy',true,'lump_factor',3); disp(nnz(p)/n)
n=1000; p = RandPartition([2 3 5 7 11 13], 100*ones(n,1),'lumpy',true,'lump_factor',4); disp(nnz(p)/n)
n=1000; p = RandPartition([2 3 5 7 11 13], 100*ones(n,1),'lumpy',true,'lump_factor',5); disp(nnz(p)/n)
n=1000; p = RandPartition([2 3 5 7 11 13], 100*ones(n,1),'lumpy',true,'lump_factor',50); disp(nnz(p)/n)
n=1000; p = RandPartition([2 3 5 7 11 13], 100*ones(n,1),'lumpy',true,'lump_factor',100); disp(nnz(p)/n)

%% 10-D Larger sample sizes (raw results hidden) for lumpiness comparisions
n=1000; p = RandPartition([2 3 5 7 11 13 ones(1,4)], 200*ones(n,1),'lumpy',false); disp(nnz(p)/n)
n=1000; p = RandPartition([2 3 5 7 11 13 ones(1,4)], 200*ones(n,1),'lumpy',true,'lump_factor',[]); disp(nnz(p)/n)
n=1000; p = RandPartition([2 3 5 7 11 13 ones(1,4)], 200*ones(n,1),'lumpy',true,'lump_factor',1); disp(nnz(p)/n)
n=1000; p = RandPartition([2 3 5 7 11 13 ones(1,4)], 200*ones(n,1),'lumpy',true,'lump_factor',2); disp(nnz(p)/n)
n=1000; p = RandPartition([2 3 5 7 11 13 ones(1,4)], 200*ones(n,1),'lumpy',true,'lump_factor',10); disp(nnz(p)/n)

%% 20-D Larger sample sizes (raw results hidden) for lumpiness comparisions
n=1000; p = RandPartition([2 3 5 7 11 13 ones(1,14)], 200*ones(n,1),'lumpy',true,'lump_factor',2); disp(nnz(p)/n)
n=1000; p = RandPartition([2 3 5 7 11 13 ones(1,14)], 200*ones(n,1),'lumpy',true,'lump_factor',20); disp(nnz(p)/n)