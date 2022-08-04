function fid = htbClose(fid)
%htbClose - Close a TEMPO HTB database file opened by htbOpen
%
% SYNOPSIS
%   err = htbClose(fid)
%
% INPUT
%   fid             fid returned by htbOpen() or htbOpenw()
%
% OUTPUT
%   err             =0 if successful, =-1 otherwise.
%
% SEE ALSO
%   htbOpenw(), htbClose(), htbGetDBCount(), htbReadHeader(), htbReadData()

err = fclose(fid);
