function mov = movieread(varargin)
% MOVIEREAD: Virtually read imphys data
% Syntax:
%   mov = movieread(tifffile)
%   mov = movieread(rawfile,headerfile)
%
% See also: IMFILE, IMPHYSLOAD.
  
% Copyright 2005 by Timothy E. Holy
  
  mov = imphysload(imfile(varargin{:}));
  