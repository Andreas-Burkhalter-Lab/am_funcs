function [ tstim ] = get_stim_time(mer,stimchan,threshold,varargin)
%GET_STIM_TIME extract stim timing from the stim channel
% This function assumes that stimulus events are indicated by superthreshold
%    values in the stim channel. 

% mer = the merec object from merecmm
% stimchan = the channel on which stimulus 
% threshold = threshold for detecting a stim event
% varargin optionally contains the desired number of stim timing to output, 
%   taken starting at the earlier superthreshold time; with no arguments,
%   only the first superthreshold time is output; input 'all' to get all
%   superthreshold times
%
% tstim = times when stim channel voltage exceeded the threshold, in scan
%    number

chanwave = mer([stimchan],[1:mer.nscans]);

if size(varargin) > 0
    if isnumeric(varargin{1})
        tstim = find(chanwave>threshold,varargin{1});
    elseif strcmp(varargin{1},'all')
        tstim = find(chanwave>threshold);
    else
        error('Last argument (number of stim times) must be numeric "all," or blank.')
    end
else
    tstim = find(chanwave>threshold,1);
end

