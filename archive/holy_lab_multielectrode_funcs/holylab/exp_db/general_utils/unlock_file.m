function suc=unlock_file(filename)
% unlock a file
% SYNTAX:
%   suc=unlock_file(filename)
   
   cmd=['dotlockfile -u ' filename '.lock'];
   [st, tt]=system(cmd);
   suc=st==0;
   