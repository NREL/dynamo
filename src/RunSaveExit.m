function RunSaveExit(command, file, run_code)
% RUNSAVEEXIT Run the given command, save workspace to file, and exit
%
% Usage:
%   RunSaveExit(command, file)
%   RunSaveExit(command, file, run_code)
%
% Warning: Exits MATLAB when done
%
% The name pretty much sums it up. This simple function runs the given
% COMMAND (passed a string), save the entire workspace to FILENAME.mat and
% exits MATLAB.
%
% This behavior is useful for long program runs, batch execution, and for
% queueing on cluster/servers
%
% Originally by Bryan Palmintier, 2011


% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-04-20 23:47  BryanP      Original version
%   2  2011-11-12 11:40  BryanP      Use v7.3 format (R2006a +) so can store >2GB if needed
%   3  2011-12-20 14:05  BryanP      Store log output from command evaluation using evalc
%   4  2012-01-25 13:30  BryanP      Use diary + temp file to capture text output without stopping realtime log
%   5  2012-01-26 11:10  BryanP      Graceful error recovery with "try". Added run_code

if nargin <3 || isempty(run_code)
    run_code = '';
else
    run_code = ['_' run_code];
end

% Save command window "diary" to a temporary file that we will later read
% in to save with the results
diary_dir = [filesep 'scratch' filesep getenv('USER') filesep 'ml_tmp_' getenv('PBS_JOBID') filesep];

% Make sure this directory exists, if not make it
if not(isdir(diary_dir))
    mkdir(diary_dir);
end

diary_file = [diary_dir 'matlab_diary.txt'];
%actually start saving to a file
diary(diary_file)

% Here we use eval() since we would like to see text in and out in real
% time for the log and diary. An alternative would be to use evalc which
% capture the output to a string, but suppresses realtime writting to
% log/diary files
%
% Enclosing in the try block allows us to save our command and diary even
% if there is a fatal error
try
    eval(command);
catch exception
    %Make the error sound
    beep
    %Print out the error statement
    disp(exception.getReport)
end

%close the diary so everything is written to it
diary('off')

% Now read in our temporary diary file so we can save this output in our
% results *.mat file
fid = fopen(diary_file);
% read the entire file as characters transpose so we have a row vector
log_eval = fread(fid, '*char')'; %#ok<NASGU> OK b/c we do use it in the eval

fclose(fid);
delete(diary_file)

%Append our
eval(['log_eval' run_code '= log_eval;']);
eval(['command' run_code '= command;']);

%Clean up
clear diary_file diary_path
if not(isempty(run_code))
    clear log_eval command
end

% Save the entire workspace
% Use v7.3 format (R2006a +) so can store >2GB if needed
save(file, '-v7.3')
exit
