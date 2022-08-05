function ip = imphys
% IMPHYS: a structure for containing imaging physiology data
% NOTE: this has been superceded by stackmm.
%
% The structure has the following fields:
%   filename: the file that an image, or it's ancestor, originated in
%   stacknum: the frame/stack number in that file
%   stacktime: the time (in seconds) at which that frame/stack was
%     acquired
%   stimulus: the valve number
%   image: the image (or vimage) itself
%
% There is a large collection of functions to manipulate these structures.
%
% See also: IMPHYSLOAD, IMPHYSFETCH, MOVIEWRITE.
  
% Copyright 2005 by Timothy E. Holy
  
  ip = struct('filename','','stacknum',0,'stacktime',[],...
              'stimulus',[],'image',[]);
  