function [sigma,x] = lightsheet_traceconfig(opparams)
% LIGHTSHEET_TRACECONFIG: do raytracing through lightsheet illumination
% Syntax:
%   [sigma,x] = lightsheet_traceconfig(opparams)
% where
%   opparams are the optical parameters, see parameters in LIGHTSHEET_PCX, plus
%     the following fields:
%       rayrgb: if set, the rays will be plotted in the current axis in
%         the color specified (3-vector [R G B]);
%       plotcomponents: if true, will plot the components in the current
%        axis;
% and
%   x is a set of positions, measured in distance from the source,
%     surrounding the ray waist;
%   sigma is the set of widths (std. dev.) of the beam as a function of
%     x.
%
% See also: LIGHTSHEET_PCX.
  
% Build the components, using plano-convex lenses
[c,r] = lightsheet_pcx(opparams);

if (isfield(opparams,'plotcomponents') & opparams.plotcomponents)
  cla
  for i = 1:length(c); ctmp = c{i}; feval(ctmp{:}); end
end

% Trace the rays, perhaps with plotting
for i = 1:length(r)
  if isfield(opparams,'rayrgb')
    rf(i) = raytrace(r(i),c,opparams.rayrgb);
  else
    rf(i) = raytrace(r(i),c);
  end
end
drawnow

%  Intensity-weighted width vs. x
[x,sigma] = raywaist(rf,[-opparams.fieldradius opparams.fieldradius]);