function [r,Lopt] = opt3d_raytrace_spherical_onaxis(r0,surfaces, ...
						  materials,drawflag)
% OPT3D_RAYTRACE_SPHERICAL_ONAXIS: 3d-raytracing for on-axis spherical surfaces
% Syntax:
%   [r,Lopt] = opt3d_raytrace_spherical_onaxis(r0,surfaces,materials,drawflag)
% where
%   r0 is a structure array of rays
%   surfaces is a structure array of surfaces, with the fields
%     c: the inverse radius of the surface (positive indicates surface
%       curves towards the +z direction, negative towards the -z direction)
%     t: the axial separation between the vertex of this surface from the
%       previous surface's vertex
%   materials is a cell array of strings describing the material to the
%     left of each surface
%   drawflag is a logical that should be true if you want to plot the
%     rays on the screen.  Optionally, it can be a vector, 1 for each
%     ray, controlling each ray individually.
% and
%   r is an nrays-by-nsurfaces+1 structure array, describing the ray at
%     each stage of the system
%   Lopt is an nrays-by-nsurfaces matrix, giving the optical path length
%     traversed at each stage of the tracing (the first column is the
%     distance to get to the first surface). The total optical path length,
%     therefore, is the sum along rows.
% 
% This is based on equations in Chapter 10 of Warren J. Smith's "Modern
% Optical Engineering"
%
% See also: PLOT_OPT3D_SURFACES.

% Copyright 2007 by Timothy E. Holy

  n_rays = length(r0);
  if isscalar(drawflag)
    drawflag = repmat(drawflag,1,n_rays);
  end
  n_surfaces = length(surfaces);
  if (surfaces(1).t ~= 0)
    error('First surface t must be 0');
  end
  % Establish a lookup table for refractive indices, since calculating them
  % is slow
  [ulambda,tmp,ulambda_lookup] = unique([r0.w]);
  [umaterials,tmp,umaterials_lookup] = unique(materials);
  nmtrx = zeros(length(umaterials),length(ulambda));
  if iscell(materials)
    for i = 1:length(umaterials)
      nmtrx(i,:) = opt_refrindx(umaterials{i},ulambda);
    end
  else
    for i = 1:length(umaterials)
      nmtrx(i,:) = umaterials(i);
    end
  end
  tcum = cumsum([surfaces.t]);
  
  % Trace the rays
  Lopt = zeros(n_rays,n_surfaces);
  r = repmat(r0(:),1,n_surfaces+1);
  for rayIndex = 1:n_rays
    % Look up the appropriate refractive indices
    n = nmtrx(umaterials_lookup,ulambda_lookup(rayIndex));
    rnew = r(rayIndex,1);
    cum_z = 0;  % Cumulative z distance travelled
    for surfIndex = 1:n_surfaces
      S = surfaces(surfIndex);
      rc = rnew;
      e = S.t * rc.e(3) - rc.x0'*rc.e;
      M1z = rc.x0(3) + e*rc.e(3) - S.t;
      M12 = rc.x0'*rc.x0 - e^2 + S.t^2 - 2*S.t*rc.x0(3);
      E1 = sqrt(rc.e(3)^2 - S.c*(S.c*M12 - 2*M1z));
      if (E1 ~= real(E1))
        if (surfIndex == 1)
          error('Ray does not strike the first surface. Did you remember to supply c as 1/R?');
        end
        % Total internal reflection, uh-oh!
        rnew.valid = false;
      end
      L = e + (S.c * M12 - 2*M1z)/(rc.e(3) + E1);
      rnew.x0 = rc.x0 + L*rc.e - [0;0;S.t];
      nratio = n(surfIndex)/n(surfIndex+1);
      Ep = sqrt(1-nratio^2*(1-E1^2));
      g1 = Ep - nratio * E1;
      rnew.e = nratio*rc.e - g1*S.c*rnew.x0 + [0;0;g1];
      % Update with respect to world-coordinates
      cum_z = cum_z + S.t;
      raynew_wc = rnew;  % The ray in world-coordinates
      raynew_wc.x0 = raynew_wc.x0 + [0;0;cum_z];
      r(rayIndex,surfIndex+1) = raynew_wc;
      Lopt(rayIndex,surfIndex+1) = n(surfIndex)*L;
      % Plot the ray, if desired
      if (drawflag(rayIndex))
        XX = [r(rayIndex,surfIndex:surfIndex+1).x0]';
        %line(XX(:,1),XX(:,2),XX(:,3),'Color',rc.rgb);
        line(XX(:,3),XX(:,2),'Color',rc.rgb);
      end
      % Plot the surface, if desired
      if any(drawflag)
        thisr = r(:,surfIndex+1);
        thisr = thisr([thisr.valid])';
        x0all = [thisr.x0];
        ymax = max(abs(x0all(2,:)));
        y = linspace(-ymax,ymax,101);
        z = S.c*y.^2./(1 + sqrt(1-S.c^2*y.^2));
        line(tcum(surfIndex)+z,y,'Color','k')
      end
    end
  end
  