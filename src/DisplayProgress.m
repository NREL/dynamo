function DisplayProgress(display_interval, n, extra_text, file_id)
% DISPLAYPROGRESS Display a series of dots to show progress
%
% DisplayProgress(display_interval, n, []) displays a series of dots. One
%     for each time n reachs an even multiple of display_interval. After 50
%     dots (or pmore precisely when n = a multiple of 50*display_interval),
%     the current value of n is displayed and the series of dots continues
%     on a new line.
%
% DisplayProgress(false, ... ) Supresses output
%
% DisplayProgress(display_interval, n, extra_text) Also displays extra_text
%     every 50 dots. extra_text is run through disp() so any thing goes;
%
% DisplayProgress(display_interval, n, [], file_id) Writes output to
% specified file handle (from fopen) rather than printing to the screen
%
% DisplayProgress(display_interval, [], ...) Uses an internal (persistent)
% counter to track the number of calls. This is particularly useful with
% parfor where a sequence may be called in any order.
%
% DisplayProgress(display_interval, [], 0) Resets the internal counter to
% equal 1.
%
% originally by Bryan Palmintier 2011

% HISTORY
% ver     date    time       who     changes made
% ---  ---------- -----  ----------- ---------------------------------------
%   1  2011-08-11 09:45  BryanP      Extracted from DP v15 & CapPlanDpOpsBatch v1
%   2  2011-12-21 15:45  BryanP      Added option to print to a file
%   3  2012-01-17 16:15  BryanP      Added internal counter
%   4  2012-01-18 20:35  BryanP      1) Return after reseting counter. 2) Fix line count
%   5  2012-02-04 02:35  BryanP      Corrected new-line issue on display/50
%   6  2012-04-23 02:35  BryanP      remove extra new line with extra_text

persistent my_n

if display_interval
    %Use our internal count if n is empty
    if nargin < 2 || isempty(n)
        % Reset count when n = [] and extra_text=0
        if isempty(my_n) || (nargin >= 3 && not(isempty(extra_text)) && extra_text == 0)
            my_n = 1;
            return
        end
        n = my_n;
        my_n = my_n + 1;
    end

    if mod(n, display_interval) == 0
        if nargin < 4 || isempty(file_id)
            %Default output to standard out (fid=1)
            file_id = 1;
        end
        fprintf(file_id,'.');

        %At end of line print current #, any extra text and return
        if mod(n, display_interval * 50) == 0
            if nargin < 3 || isempty(extra_text)
                to_print = [];
            else
                to_print = evalc('disp(extra_text)');
                to_print = strtrim(to_print);
            end
            fprintf(file_id, '%d\n%s', n, to_print);
        end
    end
end
