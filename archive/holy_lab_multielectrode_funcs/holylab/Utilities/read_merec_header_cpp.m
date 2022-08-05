function [header, headersize]=read_merec_header_cpp(fd)
% READ_MEREC_HEADER_CPP: read merec header fields.
%    This is a func similar to ReadAIHeader.
%    It reads a header created by Merec and pass a struct back to matlab
%
%    Note that: Not all fields in merec header are passed back---
%               only the fields that appear in ReadAIHeader are passed back.  
% Syntax:
%    [header, headersize]=read_merec_header_cpp(fd)
%  
% pre:
%    fd: a file handle(a.k.a file descriptor) obtained from openlfs()
%   
% post:
%    header: a matlab structure holding some fields in Merec header
%    headersize (optional) : same as header.headersize
%                    
%
% See also: OPENLFS, ReadAIHeader
