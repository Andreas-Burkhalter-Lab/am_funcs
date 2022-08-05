function WriteSonogram (fid,B,noiseThreshold,freqRange,offset,absB)
% WRITESONOGRAM: Write sonogram data to a file
%
% This function writes the sonogram DATA to a file. The data may be
% processed in blocks. The format for each block of data is the length of
% the data, all of the column numbers associated with the data, all of
% the row numbers associated with the data, all of the real parts of the
% data (complex numbers), and all of the imaginary  parts of the data
% (complex numbers).
% N.B.: the header must be written outside of this function!
%
% Syntax:
% WriteSonogram (fid,B,threshold,colOffset,absB)
% where
%    fid is the numeric file identifier
%    B is the specgram matrix
%    noiseThreshold is a value chosen to eliminate noise
%    freqRange is a 2 element vector containing the minimum and maximum
%      row numbers. The data between the min and max row numbers is
%      written to disk.
%    offset corresponds to the number of columns of B. It is incremented 
%      each block so that the data can be written in real time
%    absB is the absolute value of B (optional)
%
% See also: SOUND2SNG, READSONOGRAM.

% Copyright 2002 by Timothy E. Holy <holy@pcg.wustl.edu>
if (nargin < 6)
    absB = abs(B);
end

index = find(absB >= noiseThreshold);
[row,col] = ind2sub(size(B),index);
index2 = find(row >= freqRange(1) & row <= freqRange(2));
fwrite(fid,length(index2),'int32');
fwrite(fid,col(index2)+offset,'int32');
fwrite(fid,row(index2),'int16');
index = index(index2);
fwrite(fid,real(B(index)),'float64');
fwrite(fid,imag(B(index)),'float64');

