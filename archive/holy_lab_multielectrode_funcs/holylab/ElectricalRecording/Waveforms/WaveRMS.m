function WaveRMS
% WAVERMS: Calculate the RMS value of the waveform as a function of time
% Syntax:
%   rms = WaveRMS(filename,decfactor,trange,channellist)
% where
%   filename is a string
%   decfactor is the decimation factor
%   trange is the timerange, in seconds (default: whole file)
%   channellist is the list of channel numbers, _not_ indices
%        (default: all recorded channels)

% This is a MEX file.
