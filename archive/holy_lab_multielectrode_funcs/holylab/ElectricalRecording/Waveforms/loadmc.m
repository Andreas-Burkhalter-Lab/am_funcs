function [data,header] = loadmc(filename,time,options)
% loadmc: Load multiple-channel raw waveform data
% [data,header] = loadmc(filename)
% [data,header] = loadmc(...,timerange,options)
% Inputs:
%   filename: a string
%   fid: a file identifier, from FOPEN.
%   timerange: an optional 2-vector, indicating
%     the range [,) to load in (in seconds),
%     with 0 being the first point in the file.
%     If timerange is absent, the whole file is loaded.
%   options (optional): a structure with the following fields:
%      machfmt: machine format ('b' = big endian, etc. See FOPEN.)
%      tovolts: if true, converts to volts
%      noshow: if true, header information is not printed
% Outputs:
%        data: a matrix containing the raw waveform data (in volts)
%        Each row of the matrix corresponds to one channel
%        header: a structure containing the header information
if (nargin < 3 | ~isfield(options,'noshow'))
  options.noshow = 0;
end

[header,fid] = readheader(filename,options);
% End of header
if ~(options.noshow)
  ShowAIHeader(header);
end
fseek(fid,header.headersize,'bof');
% Bypass call to ReadBinaryData to avoid int16/uint16 issue
% in properly reading both current and legacy headers
datatype = 'uint16';
datasize = 2;
if (isfield(header,'minSample') & header.minSample < 0)
  datatype = 'int16';
end
if (isfield(header,'maxSample'))
  datasize = ceil(log2(header.maxSample)/8);
  if (datasize ~= 2)
    error('Can''t yet accomodate larger than int16s');
  end
end
nsamp = Inf;
if (nargin > 1 & ~isempty(time))
  range = round(time*header.scanrate);
  status = fseek(fid,header.numch*range(1)*datasize,'cof');
  if status
    error(ferror(fid))
  end
  nsamp = diff(range);
end
data = fread(fid,[header.numch,nsamp],datatype);
fclose(fid);
% Convert to volts
if (isfield(options,'tovolts') & options.tovolts)
  data = data*header.scalemult + header.scaleoff;
end
