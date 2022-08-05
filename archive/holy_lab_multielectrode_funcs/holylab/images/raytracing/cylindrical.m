function [t,ro] = cylindrical(r,xc,Rvec,cylax,mat1,mat2)
% The axis of symmetry in the cylinder is parallel to cylax
  P = perpproj(cylax);    % Operator to project perp. to cyl. ax.
  dx = r.x0 - xc;
  dx = P*dx(:);
  e = P*r.e(:);   % Note e will no longer be a unit vector
  % Now the rest is similar to spherical
  dx2 = sum(dx.*dx);
  e_dx = sum(e.*dx);
  R2 = sum(Rvec.*Rvec);
  % Solve quadratic equation for intersection
  sqrtarg = e_dx - dx2 + R2;
  if (sqrtarg < 0)
    % No intersection
    t = NaN;
    ro = struct;
    return
  end
  t = (-e_dx + sgn(sum(e .* Rvec)) * sqrt(sqrtarg))/sum(e.*e);
  % Now create the new ray
  xp = r.x0 + t*r.e;
  normal = P*(xp-xc);
  normal = normal/sqrt(sum(normal.*normal));
  ro = transmitted_ray(r,t,normal,mat1,mat2);
    