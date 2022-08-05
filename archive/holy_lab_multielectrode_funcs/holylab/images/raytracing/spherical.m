function [t,ro] = spherical(r,xc,Rvec,mat1,mat2)
% Rvec is a 3-vector whose amplitude specifies the radius.  It's
% direction breaks the degeneracy as to which of two intersection points
% is actually used; it points from the sphere center to the center of the
% face which is present.
% mat1 is the material on one side, mat2 is the material on the other. mat1 is the input side if the dot product (normal,r.e) is positive.
  dx = r.x0 - xc;
  dx2 = sum(dx.*dx);
  e_dx = sum(r.e.*dx);
  R2 = sum(Rvec.*Rvec);
  % Solve quadratic equation for intersection
  sqrtarg = e_dx - dx2 + R2;
  if (sqrtarg < 0)
    % No intersection
    t = NaN;
    ro = struct;
    return
  end
  t = -e_dx + sgn(sum(r.e .* Rvec)) * sqrt(sqrtarg);
  % Now create the new ray
  xp = r.x0 + t*r.e;
  normal = xp-xc;
  normal = normal/sqrt(sum(normal.*normal));
  ro = transmitted_ray(r,t,normal,mat1,mat2);
  