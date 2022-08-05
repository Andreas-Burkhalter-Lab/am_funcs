function flag = supportlegacyppnumber
% SUPPORTLEGACYPPNUMBER: a flag to determine support of old ephysplot syntax
%
% Making this a function allows the user to locally override this
% choice without changing the default setting for all users, simply by
% putting a function of the same name earlier in his/her own MATLAB
% search path.

  flag = 0;

