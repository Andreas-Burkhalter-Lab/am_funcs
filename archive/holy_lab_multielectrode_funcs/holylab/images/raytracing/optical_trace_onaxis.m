function rb = optical_trace_onaxis(rb,component)
% optical_trace_onaxis: ray tracing for layouts with a single optic axis
%
% Syntax:
%   rbOut = optical_trace_onaxis(rb,component)
% where
%   rb is the initial ray bundle (see opticalRayBundle)
%   component is a cell array of size 1-by-n, where each element is either
%     a surface (e.g., opticalQuadratic) or is a structure specifying
%     propagation, which contains the following fields:
%       t: the thickness of the material, measured along the optic axis
%         (i.e., vertex-to-vertex)
%       material: a string, e.g., 'bk7', containing the name of the optical
%         material.
% and
%   rbout is a 1-by-n+1 array of ray bundles, specifying the configuration
%   of the ray bundle at each step of tracing.
%
% See also: optical_plot_rays_2d, optical_plot_rays_3d,
% optical_plot_components.

% Copyright 2010 by Timothy E. Holy

  %% Parse arguments
  if ~iscell(component)
    error('component must be a cell array');
  end
  n_components = length(component);
  if ~isequal(rb.q,[0 0 0 1]')
    error('Optic axis must be coincident with the z-axis');
  end

  %% Establish a lookup table for refractive indices (speed optimization)
  % Collect all the different materials
  materialFlag = cellfun(@isstruct,component);
  material_structures = [component{materialFlag}];
  materials = {material_structures.material};
  % Find the unique combinations
  [ulambda,~,ulambda_lookup] = unique(rb.wavelength);
  [umaterials,~,umaterials_lookup] = unique(materials);
  nmtrx = zeros(length(umaterials),length(ulambda));
  for i = 1:length(umaterials)
    thismaterial = umaterials{i};
    if ischar(thismaterial)
      % material was supplied as the name of a type of class
      nmtrx(i,:) = opt_refrindx(thismaterial,ulambda);
    else
      % material was supplied as a number, i.e., constant refractive index
      % for all wavelengths
      nmtrx(i,:) = thismaterial;
    end
  end
  % Store the lookup index with the material structure
  materialIndex = find(materialFlag);
  for i = 1:length(materialIndex)
    component{materialIndex(i)}.material_lookup = umaterials_lookup(i);
  end

  %% Trace the rays
  rb = repmat(rb,1,n_components+1);
  for i = 1:n_components
    thiscomponent = component{i};
    if materialFlag(i)
      % Propagation
      rb(i+1) = rb(i).shift_origin(thiscomponent.t);
      thisn = nmtrx(thiscomponent.material_lookup,ulambda_lookup);
      rb(i+1) = rb(i+1).propagate_to_optic_plane(thisn);
    else
      % Refraction at a component surface
      % Find previous and next materials
      prevIndex = max(materialIndex(materialIndex < i));
      nextIndex = min(materialIndex(materialIndex > i));
      prevn = nmtrx(component{prevIndex}.material_lookup,ulambda_lookup);
      if ~isempty(nextIndex)
        nextn = nmtrx(component{nextIndex}.material_lookup,ulambda_lookup);
      else
        nextn = prevn;
      end
      % Trace through the surface
      [rb(i),rb(i+1)] = thiscomponent.trace(rb(i),prevn,nextn);
    end
  end
  
