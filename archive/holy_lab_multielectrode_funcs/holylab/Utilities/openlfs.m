function [fd, msg]=openlfs(filename)
% OPENLFS: open a file with LFS mode
% Syntax:
%    [fd, msg]=openlfs(filename)
% pre:
%    filename: the name of the file to open
% post:
%    fd: a file handle(a.k.a file descriptor, an integer) 
%    msg (optional) : system message after opening a file. Useful when
%                     error occurs. Same as what printed by perror() in c/c++.
%
% See also: CLOSELFS
