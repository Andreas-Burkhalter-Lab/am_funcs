function c = optical_layout_onaxis(varargin)
% optical_layout_onaxis: building lens systems with a single optic axis
%
% Syntax:
%   c = optical_layout_onaxis(matstruct1,surface1,matstruct2,surface2,...)
% or
%   c = optical_layout_onaxis(surface1,matstruct1,surface2,...)
% where the elements are listed in alternating orders of "materials" and
% "surfaces." The syntax for a material is a structure with fields:
%       t: the thickness of the material, measured along the optic axis
%         (i.e., vertex-to-vertex)
%       material: a string, e.g., 'bk7', containing the name of the optical
%         material.
% and surfaces are objects like opticalQuadratic.
%
% On output, c is a cell array with each component placed in the
% appropriate location, alternating with strings indicating the material
% name. The first surface is placed at [0 0 0] with an optic axis of
% [0 0 1]. See optical_shift_components and related functions if you want
% to change these parameters.
%
% You can assemble sub-assemblies by contatenation. Combining this with
% component shifting, rotation, etc, allows one to build complex, off-axis
% lens assemblies.
%
% See also: optical_plot_rays_2d, optical_plot_rays_3d,
% optical_plot_components, optical_shift_components.

% Copyright 2010 by Timothy E. Holy

  tcum = 0;
  n_components = length(varargin);
  c = {};
  for i = 1:n_components
    this_component = varargin{i};
    if isa(this_component,'opticalAperture') || iscell(this_component)
      if iscell(this_component)
        o = this_component{1}.origin;
      else
        o = this_component.origin;
      end
      this_component = optical_shift_component(this_component, [0 0 tcum]'- o);
      if iscell(this_component)
        c = [c this_component];
        tcum = tcum + this_component{end}.origin(3) - this_component{1}.origin(3);
      else
        c{end+1} = this_component;
      end
    else
      c{end+1} = this_component.material;
      tcum = tcum + this_component.t;
    end
  end
