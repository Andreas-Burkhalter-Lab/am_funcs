function [onsandoffs,numofstims] = onsandoffs(basename)
%onsandoffs is a parsing function that takes the header file for an Imagine
%experiment, does indexing gymnastics, and returns an array of cells
%containing matrices that indicate the on and off times for each stimulus
%pulse. The cell array is indexed by the valve number of the stimulus and
%the matrices each have two columns.  The left column indicates the very
%last flush stack before the stimulus comes on (which we choose as the 
%indicator for the valve 'on moment' because the valve turns on immediately 
%after this stack is completed) and the right column indicates the last 
%'valve on' stack. The rows of each matrix index trials for a given stimulus. 
%
%onsandoffs is the returned cell array of matrices with indexed transitions
%numofstims is just a scalar indicator of the total number of stimuli presented
%
%There, that should be clear as mud.
%
%syntax:
%       [onsandoffs,numofstims] = onsandoffs(basename)
%where:
%      basename is the original name of the imagine experiment
%      onsandoffs is a cell array with dimension described above
%      numofstims tells you the total number of stimulus presentations
%
%Copywrite 2006 by Terrence Holekamp

h = imreadheader(basename);
onsandoffs = cell(length(h.stim_labels),1);
stimnum = nan(length(h.stim_labels),1);
for i = 1:length(h.stim_labels)
    stackind = [0 find(h.stim_lookup == i)']';
    stackind(stackind > h.nstacks) = [];
    transindex = diff(stackind,1,1);
    transindex(transindex == 1) = 0;
    cliffs = find(transindex);
    onindex = (stackind(cliffs+1))-1;
    offindex = (stackind(cliffs));
    offindex(1) = [];
    offindex = [offindex' stackind(end)']';
    onsandoffs{i} = [onindex offindex];
    stimnum(i) = length(onindex);
end

numofstims = nansum(stimnum);

















