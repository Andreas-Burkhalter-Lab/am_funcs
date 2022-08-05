 function [sng,header] = ReadSonogram(file,timeRange)
% READSONOGRAM: Read sonogram data from a file
%
% This function reads the header and the sonogram data from a file. It can 
% read the entire file or a specified portion of the file indicated by a 
% time range. The read may be processed in blocks.  The format for each block
% of data read is the length of the data, all of the column numbers associated 
% with the data, all of the row numbers associated with the data, all of the 
% real parts of the data (complex numbers), and all of the imaginary parts of the
% data (complex numbers).
%
% Syntax:
% [sng,header] = ReadSonogram(file,timeRange)
% where
%    file is either a string and is treated as a filename or
%      is numeric and is treated as a file identifier
%    timeRange is vector containing the start and stop times (which correspond
%      to the aquisition time) of the read (optional). If a timeRange is not 
%      provided, the entire file is read.
% and
%    sng is a sparse matrix; each entry contains the row number, the column number,
%      the real part of the complex number, and the imaginary part of the 
%      complex number
%    header is a structure that returns the information in the header
%
% See also: WRITESONOGRAM

% Copyright 2002 by Timothy E. Holy <holy@pcg.wustl.edu>
  if (ischar(file))
    [fid,message] = fopen(file,'r');
    if (fid < 1)
      error(message);
    end
  elseif (isnumeric(file))
    fid = file;
  else
    error(['Do not recognize input ',file]);
  end
  [header,headersize] = ReadAIHeaderWU1(fid);
  % If timeRange isn't provided (nargin < 2), then
  % default to entire timerange
  if (nargin < 2)
    timeRange = [0 header.tacq];
  end
  if (timeRange(2) > header.tacq)
    errordlg('Time range exceeds total aquistion time', 'Input Error');
    return;
  end
  % Convert to column numbers from the times
  nColsPerSec = (header.columnTotal/header.nblocks)/(header.tacq/header.nblocks);
  colStartRead = nColsPerSec * timeRange(1) + 1;
  colStopRead = nColsPerSec * timeRange(2) + 1;
  
  nbufs = 0; 
  lastColumn = 0;
  colAccum = cell(1,header.nblocks);
  rowAccum = colAccum;
  realAccum = colAccum;
  imagAccum = colAccum;
  numberOfBytes = 2+2*8;
  
  while (nbufs < header.nblocks && lastColumn < colStopRead)
    length = fread(fid,1,'int32');
    col = fread(fid,length,'int32');
    index = find(col >= colStartRead & col <= colStopRead);
    if ~isempty(index)
      colAccum{nbufs+1} = col(index);
      row = fread(fid,length,'int16');
      rowAccum{nbufs+1} = row(index);
      realB = fread(fid,length,'float64');
      realAccum{nbufs+1} = realB(index);
      imagB = fread(fid,length,'float64'); 
      imagAccum{nbufs+1} = imagB(index);
    else
      fseek(fid,length*numberOfBytes,0);
    end
    if (length == 0)
      lastColumn = 0;
    else   
      lastColumn = col(end);
    end
    nbufs = nbufs + 1;
  end
  colAccum = cat(1,colAccum{:});
  rowAccum = cat(1,rowAccum{:});
  realAccum = cat(1,realAccum{:});
  imagAccum = cat(1,imagAccum{:});

  Belements = realAccum + imagAccum*1i;

  % Convert to sparse matrix
  sng = sparse(rowAccum,colAccum,Belements,header.nfreq,header.columnTotal);
  
  % Close file if "file" is a string
  if (ischar(file))
    status = fclose(fid);
    if (status < 0)
      error('File did not close');
    end
  end

