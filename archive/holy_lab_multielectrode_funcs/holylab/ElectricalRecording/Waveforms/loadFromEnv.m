function [data,header] = loadFromEnv(filename,time)
% [data,header] = loadFromEnv(filename,timerange)
% Loads envelope data (.env files)
% data contains a matrix with 2*numCh rows,
% row 1/2 are min/max of first channel,
% row 3/4 are min/max of second channel, etc.
% See loadenv to get envelopes from .bin files
% (or better, run Parse)
% See fillEnv for plotting

%if ischar(filename)
%  [fid,message] = fopen(filename,'r');
%  if (fid < 1)
%    error(message)
%  end
%elseif isnumeric(filename)
%  fid = filename;
%else
%  error('Invalid filename/fid input');
%end
[header,fid] = readheader(filename);
% End of header
ShowAIHeader(header)
fprintf('Decimation factor: %f\n',header.decimate)
if (nargin > 1)
  range = round(time*header.scanrate/header.decimate);
else
  range = [0,header.nscans];
end
data = ReadBinaryData(fid,2*header.numch,range);
if ~isnumeric(filename)
  fclose(fid);
end
% Convert to volts
data = data*header.scalemult + header.scaleoff;
%for i = 1:header.numch
%        data(2*i-1,:) =data(2*i-1,:)*header.scalemult + header.scaleoff;
%        data(2*i,:) = data(2*i,:)*header.scalemult + header.scaleoff;
%end
