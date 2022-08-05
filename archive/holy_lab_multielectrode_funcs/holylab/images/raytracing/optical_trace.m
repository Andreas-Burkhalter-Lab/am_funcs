function [rb,varargout] = optical_trace(rb,component,options)
% optical_trace: ray tracing through component layouts
%
% Syntax:
%   rbOut = optical_trace(rb,component)
%   [rbOut,extraOutputs] = optical_trace(rb,component,options)
% On output,
%   rbOut is a 1-by-n array of ray bundles (see opticalRayBundle),
%   specifying the configuration of the ray bundle at each surface
% These are the inputs:
%   rb is the initial ray bundle 
%   component is a cell array of size 1-by-n, containing alternating
%     materials and surfaces. Each component must be placed appropriately,
%     which is done with other tools (see optical_layout_onaxis, for
%     example). The component list must start with a material name.
%   options is a structure which may have the following fields:
%     tofocus: if present, continues tracing the rays beyond the last
%       surface to a point of focus. The value of the field 'tofocus'
%       should be
%              options.tofocus = normal(:);
%       where 'normal' is a vector specifying the normal vector of the
%       plane in which you wish to achieve focus.  If you want to define
%       focus as being "narrowest along just the axis 'v'" rather than over
%       the whole focal plane, then set
%             options.tofocus = [normal(:) v(:)];
%       When 'tofocus' is supplied, the first two "extraOutputs" are
%       mean_coef and cov_coef, based from the plane of focus (see
%       opticalRayBundle.mean_covariance).
%
% See also: opticalRayBundle, optical_layout_onaxis, optical_plot_rays_2d,
% optical_plot_rays_3d, optical_plot_components_2d, optical_plot_components_3d.

% Copyright 2010 by Timothy E. Holy

  %% Parse arguments
  if ~iscell(component)
    error('component must be a cell array');
  end
  if (nargin < 3)
    options = struct;
  end

  %% Establish a lookup table for refractive indices (speed optimization)
  % Collect all the different materials
  materialFlag = cellfun(@ischar,component) | cellfun(@isnumeric,component);
  if ~materialFlag(1)
    error('First item in component list must be a material');
  end
  if any(diff(materialFlag) == 0)
    error('component must be an alternating list of materials and surfaces');
  end
  materials = component(materialFlag);
  if (materialFlag(end) == false)
    % If no final materials is specified, continue with the last material
    % (no refraction at final surface)
    materials{end+1} = materials{end};
  end
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
%   % Store the lookup index with the material structure
%   materialIndex = find(materialFlag);
%   for i = 1:length(materialIndex)
%     component{materialIndex(i)}.material_lookup = umaterials_lookup(i);
%   end

  %% Trace the rays
  surfaceIndex = find(~materialFlag);
  n_surfaces = length(surfaceIndex);
  rb = repmat(rb,1,n_surfaces+1);
  for i = 1:length(surfaceIndex)
    this_component = component{surfaceIndex(i)};
    npre = nmtrx(umaterials_lookup(i),ulambda_lookup);
    npost = nmtrx(umaterials_lookup(i+1),ulambda_lookup);
    rb(i+1) = this_component.trace(rb(i),npre,npost);
  end

  %% Trace to focus, if desired
  if isfield(options,'tofocus')
    if (numel(options.tofocus) == 3)
      options.tofocus = options.tofocus(:);
    end
    normal = options.tofocus(:,1);
    [mean_coef,cov_coef] = rb(end).mean_covariance(normal);
    c = [0 0 0];
    if (size(options.tofocus,2) == 1)
      % Full focus
      for i = 1:3
        c(i) = trace(cov_coef(:,:,i));
      end
    else
      % Focus along axis 'v'
      v = options.tofocus(:,2);
      for i = 1:3
        c(i) = v'*cov_coef(:,:,i)*v;
      end
    end
    if (c(1) == 0)
      error('Beam is collimated, there is no focus');
    end
    % Find the position of the focus
    zmin = -c(2)/(2*c(1));
    pf = zmin*mean_coef(:,1) + mean_coef(:,2);
    % Propagate to the plane containing pf
    screen = opticalAperture(Inf,pf);
    rb(end+1) = screen.trace(rb(end),npost);
    % Return mean and covariance, re-centered on the focal plane
    mean_coef(:,2) = pf;
    cov_coef(:,:,3) = zmin^2*cov_coef(:,:,1) + zmin*cov_coef(:,:,2) + cov_coef(:,:,3);
    cov_coef(:,:,2) = 2*zmin*cov_coef(:,:,1) + cov_coef(:,:,2);
    varargout = {mean_coef,cov_coef};
  end
