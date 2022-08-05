function wave=readuint16lfs(fd, nChannels, scanRange, datastart)
% READINT16LFS: read a block of uint16's into wave
% Syntax:
%    wave=readuint16lfs(fd, nChannels, scanRange, datastart)
%  
% pre:
%    fd: a file handle(a.k.a file descriptor, an integer) obtained by openlfs()
%    nChannels: #channels in each scan
%    scanRange: a 2-element vector [a,b] defines the range to read. Note
%               that the range includes both a and b.
%               Also a, b takes values in [0, nTotalScans).
%    datastart: a number tells where Merec data begins, i.e. the Merec
%               header size.  
%  
% post:
%    wave: a nChannels X (scanRange(2)-scanRange(1)+1) matrix whose elts
%          are int16's 
%
% examples:
%    see snippetfile.m 
%  
% See also: OPENLFS
% this is the old slow code:   
% wave=uint16_c(readint16lfs(fd, nChannels, scanRange, datastart) );
