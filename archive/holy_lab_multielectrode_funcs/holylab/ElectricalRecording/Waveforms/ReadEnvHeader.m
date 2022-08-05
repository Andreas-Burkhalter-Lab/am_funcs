function [header,headersize] = ReadEnvHeader(fid)
% [header,headersize] = ReadEnvHeader(fid)
% Read the header from .env files
[header,headersize] = ReadAIHeaderHarvard(fid);
if (header.type ~= 4)
        error('File is not .env type');
end
header.versionEnv = fread(fid,1,'int16');
header.decimate = fread(fid,1,'int32');
header.nscans = header.nscans*header.decimate;  % TEH 05-08-2003
