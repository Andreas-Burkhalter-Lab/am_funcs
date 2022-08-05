function [header,headersize] = ReadAIHeaderWU1(file)
% READAIHEADERWU1: Read the header of raw waveform files (AI = analog input)
% 
% This function reads the header of many different file types. It reads 
% version 3 of the universal AI header type and version 1 of the other 
% header types.  A magic number is also read at the beginning of each 
% header to ensure proper endian status on a particular machine and to
% check if the header should be read by a legacy function.
%
% Syntax:
% [header,headersize] = ReadAIHeaderWU1(file)
% where
%    file is either a string and is treated as a filename or
%      is numeric and is treated as a file identifier
% and
%    header is a structure that returns the information in the header
%    headersize is the size of the header in bytes
%
% See also: WRITEAIHEADERWU1

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
magicNum = fread(fid,1,'int32');
if (magicNum == 992732356)
  headersize = fread(fid,1,'uint32');
  header.headersize = headersize;
  type = fread(fid,1,'int32');
  version = fread(fid,1,'int16');
  if (bitget(type,1) & version == 3)
    % Read universal AI header type
    header.nscans = fread(fid,1,'uint32');
    numChannels = fread(fid,1,'int32');
    header.numch = numChannels;
    header.channels = fread(fid,numChannels,'int16');
    if (size(header.channels,2) == 1)
      header.channels = header.channels';
    end
    header.scanrate = fread(fid,1,'float32');
    header.scalemult = fread(fid,1,'float32');
    header.scaleoff = fread(fid,1,'float32');
    header.voltageMin = fread(fid,1,'int16');
    header.voltageMax = fread(fid,1,'int16');
    header.date = ReadLVString(fid);
    header.time = ReadLVString(fid);
    header.usrhdr = ReadLVString(fid);
  end
  if (bitget(type,4))
    % Read sonogram header type 
    version = fread(fid,1,'int16');
    header.nfreq = fread(fid,1,'int32');
    header.columnTotal = fread(fid,1,'int32');
    header.threshold = fread(fid,1,'float32');
    header.nblocks = fread(fid,1,'int32');
    header.tacq = fread(fid,1,'float32');
    header.freqMin = fread(fid,1,'int32');
    header.freqMax = fread(fid,1,'int32');
  end
  if (bitget(type,5))
    % Read proximity detection header type
    version = fread(fid,1,'int16');
    header.numTransitions = fread(fid,1,'int16');
  end
  % Set file position indicator to the end of the header
  fseek(fid,headersize,-1);
  header.flag = 'WU1';
  
elseif (magicNum == -991679685)
  error(['File has non-native endian status.  Call FOPEN with the reverse' ...
      ' byte order (see FOPEN help)']);
else
  warning('Assuming legacy data file; proceding with old-style header')
  fseek(fid,0,-1);
  [header,headersize] = ReadAIHeaderHarvard(fid);
end
% Close file if "file" is a string
if (ischar(file))
  status = fclose(fid);
  if (status < 0)
    error('File did not close');
  end
end
