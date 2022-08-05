function [fid,npts,numch] = posaifile(filename,trange)
% POSAIFILE: position analog input (AI) file to right location
% [fid,npts,numch] = posaifile(filename,tstart)
% where
%   filename is a string containing the name of the file
%   trange is the timerange [tstart tend] (in seconds) of the data
%     you wish to read
% and
%   fid is the file id
%   npts is the total number of scans to read
%   numch is the number of channels recorded in the file
%[fid,message] = fopen(filename,'r');
  [fid,message] = fopen(filename,'r','b');
if (fid < 1)
        error(message)
end
header = ReadAIHeader(fid);
status = fseek(fid,header.headersize,'bof');        % Move to the beginning of data
if (status)
        error(ferror(fid))
end
if (~isempty(trange))
        irange = round(trange*header.scanrate);
else
        irange = [0 header.nscans-1];
end
npts = diff(irange)+1;
status = fseek(fid,header.numch*irange(1)*2,'cof');        % Move to the beginning of timerange
if (status)
        error(ferror(fid))
end
numch = header.numch;
