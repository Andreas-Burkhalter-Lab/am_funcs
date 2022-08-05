function e = generate_rays(objects,n_rays,lens)
% GENERATE_RAYS: creates a collection of rays that passes through a lens
% Syntax:
%   e = generate_rays(objects,n_rays,lens)
% where
%   objects is a d-by-n_objects matrix, where d is the dimensionality of
%     space (3, or 2 if you want only meridional rays/2d optics). Each
%     column gives the coordinates of a source of rays (an "object point")
%     in the object space, relative to the axial point of the focal plane.
%     (Note that objects do not have to be in the focal plane. The dth
%     coordinate represents the displacement along the optic axis, so
%     should be zero if you want to specify positions in the focal plane.
%     Positive values indicate positions closer to the objective than the
%     focal plane.)
%   n_rays is a scalar, giving a target for the number of generated rays.
%     This target may not be precisely met, because some generated rays
%     will not pass through the back aperture, and will be excluded from
%     the output
%   lens is a structure with the following fields:
%     NA: the numerical aperture of the "objective"
%     n_i: the index of refraction of the immersion medium in the object
%       space
%     f: a 2-vector containing the [front back] focal lengths of the
%       objective. The sign convention is such that front < 0, back > 0. If
%       you don't know what to specify here: fback = 180/mag (for Olympus
%       tube lenses with 180mm focal length), where mag is the stated
%       magnification of the lens. Then, use ffront = -n_i*fback (which
%       models the objective as a single-component ideal lens).
%     yrange (optional): a 2-vector, [min max], in normalized units that
%       specifies the window of allowed y-coordinates in the back focal
%       plane ([-1 1] would allow all to pass, a more restrictive range
%       would correspond to blocking a portion of the objective)
% and
%   e is a cell array, with length equal to the number of object points,
%     where each element is a d-by-n matrix. Each column gives the
%     components of a unit vector that passes from the corresponding source
%     point through the objective's back aperture. n may be different for
%     different source points. The raytracing is done assuming an ideal
%     sine-lens objective.

  [d,n_objects] = size(objects);
  a = lens.NA/lens.n_i;  % the normalized back aperture size
  r_app = abs(a*lens.f(1));  % the actual back aperture size
  if (d == 2)
    n_rays_lin = n_rays;
  elseif (d == 3)
    n_rays_lin = ceil(sqrt((4/pi)*n_rays));
  else
    error('d > 3 not implemented yet');
  end
  n_rays_created = n_rays_lin^(d-1);
  ecoords = cell(1,d-1);
  egrid = cell(1,d-1);
  e = cell(1,n_objects);
  for objIndex = 1:n_objects
    % For each object, compute the range of angles that will be probed to
    % insure that the entire back aperture could be filled
    x0tilde = objects(1:d-1,objIndex);
    norm_objpos = x0tilde/lens.f(1);
    norm_disp = sqrt(sum(norm_objpos.^2));
    emax = norm_disp + a;  % upper bound on sine of angle relative to optic axis
    % Generate the rays as an evenly-spaced grid
    ey = linspace(-1,1,n_rays_lin) * emax;
    if (d > 2)
      for i = 1:d-1
        egrid{i} = ey;
      end
      [ecoords{:}] = ndgrid(egrid{:});
    else
      ecoords{1} = ey;
    end
    eout = zeros(d,n_rays_created);
    for i = 1:d-1
      eout(i,:) = ecoords{i}(:)';
    end
    eout(d,:) = sqrt(1 - sum(eout(1:d-1,:).^2,1));
    % For each ray, check to see if it passes through the aperture
    % First we propagate to the focal plane
    x0tp = repmat(x0tilde,1,n_rays_created) - repmat(objects(d,objIndex)./eout(d,:),d-1,1) .* eout(1:d-1,:);
    % Now check whether rays make it through the back aperture
    x1t = x0tp - lens.f(1) * eout(1:d-1,:);
    valid = sum(x1t.^2,1) <= r_app^2;
    if isfield(lens,'yrange')
      valid = valid & x1t(2,:) > lens.yrange(1)*r_app & x1t(2,:) < lens.yrange(2)*r_app;
    end
    e{objIndex} = eout(:,valid);
  end
