function out = rms(list)
% RMS: Root Mean Square of a vector.
%
% Basic RMS function to work around the fact that RMS is only included in
% the signal processing toolbox.
%
% Adapted from http://rosettacode.org/wiki/Averages/Root_mean_square#MATLAB
%
% If passed a matrix, will work by column
    out = sqrt(mean(list.^2));
end