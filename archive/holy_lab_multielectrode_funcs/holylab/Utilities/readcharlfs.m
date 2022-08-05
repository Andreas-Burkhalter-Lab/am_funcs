function buf=readcharlfs(fd, unitsize, unitRange, datastart)
% READCHARLFS: read a block of chars into buf
% Syntax:
%    buf=readcharlfs(fd, unitsize, unitRange, datastart)
%  
% pre:
%    fd: a file handle(a.k.a file descriptor, an integer) obtained by openlfs()
%    unitsize: #bytes in each unit
%    unitRange: a 2-element vector [a,b] defines the range to read. Note
%               that the range includes both a and b.
%              
%    datastart: a number tells where data begins, i.e. the first byte read at:
%               datastart+firstunit*unitsize;   
%  
% post:
%    buf: a unitsize X (unitRange(2)-unitRange(1)+1) matrix whose elts
%          are chars 
% samples:
%    buf=readcharlfs(fd, 1, [0 1023], 0);   
% See also: OPENLFS
   
   buf=char(readint8lfs(fd, unitsize, unitRange, datastart));
