function T = ocpi_array2native(instrument_params, ds, modestr, restrict_schedule)
% ocpi_array2native: create tform to convert array coordinates to physical units
% Here "native" means that the xy plane of the two coordinate systems are
% parallel (i.e., camera-native).
%
% Syntax:
%   T = ocpi_array2native(instrument_params, ds, modestr)
% where
%   instrument_params is a structure of the format produced by ocpi_params;
%   ds is the displacement along the scanning axis between frames (microns)

% See also: ocpi_params.

% Copyright 2012 by Timothy E. Holy

  mode = find(strcmpi(modestr,{'2d','3d'}));
  if isempty(mode)
    error(['Mode ' modestr ' not recognized, use 2D or 3D'])
  end
  dxy = instrument_params.pixel_xy_microns;
  theta = instrument_params.angle_scan_optax;
  have_restrict_schedule = nargin > 3 && ~isempty(restrict_schedule);
  
  if have_restrict_schedule
    scale = prod(2.^restrict_schedule,1);
  end
  
  A = [dxy -ds*sin(theta); 0 ds*cos(theta)];
  icoord = abs(instrument_params.illumsource_coord);
  c = [icoord 3];
  if mode == 1
    % 2d case (i.e., a y-slice, i.e., perpendicular to light sheet horizontal plane)
    if have_restrict_schedule
      A = A * diag(scale(c));
    end
    T = maketform('affine',[A; 0 0]);
  else
    %3d case
    A3 = zeros(3,3);
    A3(c,c) = A;
    ncoord = setdiff([1 2],icoord);
    A3(ncoord,ncoord) = dxy;
    if have_restrict_schedule
      A3 = A3 * diag(scale);
    end
    T = maketform('affine',[A3; 0 0 0]);
  end
  
end
