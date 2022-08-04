function fid = htbOpenw(filename)
%htbOpenw - Open a TEMPO HTB database file for read/write access
%
% The file is opened for reading or appending to the end
%
% SYNOPSIS
%   fid = htbOpenw(filename)
%
% INPUT
%   filename        Name of HTB file to open
%
% OUTPUT
%   fid             Matlab file identifier
%                   =-1 if there was an error opening HTB file
%
% SEE ALSO
%   htbOpenw(), htbClose(), htbGetDBCount(), htbReadHeader(), htbReadData()

fid = fopen(filename, 'a+');
return;