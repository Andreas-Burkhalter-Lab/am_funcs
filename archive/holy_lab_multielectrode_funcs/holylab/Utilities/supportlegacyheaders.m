function slh = supportlegacyheaders
% SUPPORTLEGACYHEADERS: a flag determining data format compatibility
%
% Making this a function allows the user to locally override this
% choice without changing the default setting for all users, simply by
% putting a function of the same name earlier in his/her own MATLAB
% search path.

  slh = 1;   % After backward-compatibility testing is done, this will
             % default to 0.
