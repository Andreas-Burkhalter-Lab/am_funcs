function params = ocpi_params(instrument)
% ocpi_params: return geometric parameters of various OCPI microscopes at WashU
%
% Syntax:
%   params = ocpi_params(instrument)
% where
%   instrument is a string, 'OCPI1', 'HS-OCPI'
% and
%   params is a structure whose field names should be largely
%     self-explanatory. For example, 'angle_scan_optax' means the angle
%     between the scan direction and the optical axis of the objective (in
%     radians).
%     The illumsource_coord describes where the illumination comes from. It
%     is a signed integer; -2 would mean that it propagates along the 2nd
%     image coordinate, coming from the "1" side (rather than the "end"
%     side).

% Copyright 2012 by Timothy E. Holy

  params.instrument = instrument;
  switch instrument
    case 'OCPI1'
      params.numAper = 0.5;
      params.pixel_xy_microns = 0.71;
      params.angle_scan_optax = 0;
      params.angle_scan_horiz = 40*pi/180;  % beware, can be changed
      params.illumsource_coord = 2;  % from positive, 2nd coordinate
    case 'HS-OCPI'
      params.numAper = 0.5;
      params.pixel_xy_microns = 0.325;
      params.angle_scan_optax = 52.4*pi/180;
      params.angle_scan_horiz = 0;
      params.illumsource_coord = 2;
    otherwise
      error('Setting unknown')
  end
end
