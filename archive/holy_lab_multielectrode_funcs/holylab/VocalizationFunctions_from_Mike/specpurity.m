    function [totalPower, maxPower, numFreqsAccum] = specpurity(powerspec)
% SPECPURITY: Distinguishes singing from non-singing
%
% This function computes the total power, the maximum power, and the
% the number of frequencies for each time bin (column). This information
% is then used to distinguish singing from non-singing.
%
% Syntax:
% [totalPower, maxPower, numFreqsAccum] = specpurity(powerspec)
% where
%    powerspec is a sparse matrix of the sonogram data squared
% and
%    totalPower is a vector of the sums of the sonogram data in each column
%    maxPower is a vector of the maximum values of the sonogram data in each column
%    numFreqsAccum is a vector of the number of frequencies in each column

% Copyright 2002 by Timothy E. Holy <holy@pcg.wustl.edu>

totalPower = sum(powerspec);
maxPower = max(powerspec);
[row,col] = find(powerspec);
index = find(diff(col)~=0);

firstEntry = index(1);
numFreqs = diff(index);
lastEntry = length(col)-index(end);

numFreqsAccum = zeros(1,length(totalPower));
numFreqsAccum([col(index); col(end)]) = [firstEntry; numFreqs; lastEntry];
