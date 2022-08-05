function [angles,options] = rotatestack_prepare_imagine(header,lightdir,angle)
% ROTATESTACK_PREPARE_IMAGINE: set up data for stack rotation with IMAGINE files
% Syntax:
%   [angles,options] = rotatestack_prepare_imagine(header,lightdir,angle)
% where
%   header is the header structure from reading the imagine file (e.g.,
%     stackmm.header or imreadheader)
%   lightdir is a string, 'left','right','top','bottom', indicating the
%     direction from which the illumination is coming.  This should be
%     evaluated using IMAGESC, _not_ STACKVIEWER (which performs a
%     rotation on the raw data).
%   angle is an angle in radians, most commonly the tilt angle of the
%     microscope relative to horizontal.  This can be ommitted if the
%     .imagine file has an "angle" entry.
% and
%   angles is a 3-vector suitable for input to ROTATESTACK.
%   options is an options structure suitable for input to
%     ROTATESTACK. Note, however, that no cropping is specified.  See
%     ROTATESTACK for information.
%
% See also: ROTATESTACK.
  
% Copyright 2009 by Timothy E. Holy
  
  options.pixel_spacing = [[1 1]*header.um_per_pixel_xy ...
		    diff(header.piezo_start_stop)/(header.depth-1)];
  if (isfield(header,'angle') && nargin < 3)
    angle = pi*header.angle/180;
  end
  switch lower(lightdir(1))
   case 'l'
    angles = [angle 0 0];
   case 'r'
    angles = [-angle 0 0];
   case 't'
    angles = [0 angle 0];
   case 'b'
    angles = [0 -angle 0];
  end
end