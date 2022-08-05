function ro = transmitted_ray(r,normal_mtrx,mat1,mat2)
% TRANSMITTED_RAY: calculate the direction of transmitted ray
% Syntax:
%   ro = transmitted_ray(r,normal,mat1,mat2)
% where
%   r is a structure array of input rays
%   normal is the surface normal at the point of intersection, an
%     n_dims-by-n_rays matrix, or a single n_dims-by-1 vector if that
%     normal is identical for all rays
%   mat1,mat2 are materials 1 and 2. mat1 is the "input" side if the dot
%     product of normal with r.e is positive.
%
% and
%   ro is the output ray.

  ro = r;
  if isequal(mat1,mat2)
    % The materials are identical, so there is no refraction
    return;
  end
  n_dims = length(r(1).x0);
  e_all = [r.e];
  [w,tmp,wIndex] = unique([r.w]);
  nratio_all = opt_refrindx(mat1,w)./opt_refrindx(mat2,w);
  if (size(normal_mtrx,2) == 1)
    normal = normal_mtrx;
    P = perpproj(normal);
    N = diag(normal);
  end
  n_rays = length(r);
  for i = 1:n_rays
    e = e_all(:,i);
    if (size(normal_mtrx,2) > 1)
      normal = normal_mtrx(:,i);
      P = perpproj(normal);
      N = diag(normal);
    end
    nratio = nratio_all(wIndex(i));
    % Calculate the ratio of indices of refraction
    if (sum(normal.*e) < 0)
      nratio = 1/nratio;
    end
    % Calculate the direction of propagation
    e2 = nratio * (P*e);  % Snell's law
    if (norm(e2) <= 1)
      % refraction
      ro(i).e = e2 + sign(sum(e .* normal)) * sqrt(1-sum(e2.*e2))*normal;
      % Ignore intensities (i.e., reflected ray) for now
    else
      % total internal reflection
      ro(i).e = e - 2*sum(e.*normal)*normal;
    end
  end
