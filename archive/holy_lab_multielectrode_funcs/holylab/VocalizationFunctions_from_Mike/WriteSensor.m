function [numTransitions,transitions] = WriteSensor(fid,voltage,finalVoltage,threshold,offset,numTransitions)
% WRITESENSOR:  Write the proximity detection data to a file
%
% This function writes proximity detection data to a file.  The data may be
% processed in blocks. Each block of data consists of an index of low to high
% and high to low transitions sorted in ascending order. The first block of data
% also includes the initial voltage which will determine if the first transition 
% is low to high or high to low. The remaining transitions can then be deduced from
% this information. The header must be written outside of this function.
%
% Syntax:
% numTransitions = WriteSensor(fid,voltage,finalVoltage,thresholds,offset,numTransitions)
% where
%    fid is the numeric file identifier.
%    voltage is the voltage of the proximity detector calculated from from the ADD card
%    finalVoltage is equal to the initial voltage on the first block. On all other blocks,
%      the finalVoltage is equal to the the final voltage on the previous block.
%    threshold is a 2 element vector containing the low and high thresholds of the proximity detector.
%    offset corresponds to the block size. It is incremented each block so that the data can 
%      be written in real time.
%    numTransitions contains the number of transitions from low to high or high to low.  It is
%      outputed and inputed to keep a transition total for all of the blocks.
%
% See also: READSENSOR

lowThreshold = threshold(1);
highThreshold = threshold(2);

% Write initial voltage to disk
if (offset == 0)
    initialVoltage = finalVoltage;
    fwrite(fid,initialVoltage,'float64');
end

% Find low to high transitions
LtoHtransitions = find([finalVoltage voltage(1:end-1)] < lowThreshold & voltage >= highThreshold);
LtoH = LtoHtransitions + offset;
numTransitions = numTransitions + length(LtoHtransitions);

% Find high to low transitions
HtoLtransitions = find([finalVoltage voltage(1:end-1)] >= highThreshold & voltage < lowThreshold);
HtoL = HtoLtransitions + offset;
numTransitions = numTransitions + length(HtoLtransitions);

% Sort transitions in ascending order and write to disk
transitions = sort([LtoH HtoL]);
fwrite(fid,transitions,'int32');




