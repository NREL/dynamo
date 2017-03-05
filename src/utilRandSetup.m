function adp_opt = utilRandSetup(adp_opt)
% UTILRANDSETUP Setup consistent random numbers if desired
%
%   adp_opt = utilRandSetup(adp_opt)
%
% Abstracts out configuration of random numbers.
% 
% For use with the adp* family of functions.

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2017-03-04 23:00  BryanP      Adapted from adpSBI1 v23

if adp_opt.fix_rand && not(adp_opt.fix_rand_is_done)
    disp('  RAND RESET: adpSBI using fixed random number stream')
    % Reset the random number generator
    rng('default');

    %Old (pre-R2011a) syntax
	%rand_num_sequence = RandStream('mt19937ar');
	%RandStream.setDefaultStream(rand_num_sequence)

    %Indicate that we have already reset the random number generator to
    %prevent sub-functions from resetting again
    adp_opt.fix_rand_is_done = true;
end