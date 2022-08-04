function fid = htbOpen(filename)
%htbOpen - Open a TEMPO HTB database file for read only access
%
% SYNOPSIS
%   fid = htbOpen(filename)
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

fid = fopen(filename, 'r');