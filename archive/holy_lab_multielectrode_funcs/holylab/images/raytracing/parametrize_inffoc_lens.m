function [G,alpha,lens_out] = parametrize_inffoc_lens(phi,xfin,sfin,dxinf,lens)
  %% Input validation
  if (~iscell(dxinf) || ~iscell(sfin))
    error(['Inputs 2-4 must be cell arrays, with one entry per bundle of' ...
	   ' parallel rays']);
  end
  n_phi = length(phi);
  if (n_phi ~= numel(dxinf) || n_phi ~= numel(sfin))
    error('Inputs 2-4 must all have the same number of entries');
  end
%   if (~isnumeric(zN) || ~isequal(size(zN),[n_phi 2]))
%     error(['zN must be numeric, of size n-by-2, where n is the number of' ...
% 	   ' ray bundles']);
%   end
  d = size(dxinf{1},1);
  phi_zindx = find(phi == 0);
  if length(phi_zindx) ~= 1
    error('Must supply a single ray bundle with phi = 0');
  end
  if (d == 2)
    z_indx = 1;  % optic axis is x in 2d optics
  elseif (d == 3)
    z_indx = 3;  % optic axis is z in 3d optics
  else
    error('x and s must be either 2- or 3-dimensional');
  end
  nonz_indx = setdiff(1:d,z_indx);
  
  % Determine convergence locations
  for i = 1:n_phi
    [err,x0(:,i)] = ray_convergence_error(xfin{i},sfin{i});
  end
  %zF = x0(z_indx,:) - x0(z_indx,phi_zindx);
  
  % Determine G1 and alpha
  f = [lens.fL lens.fR];
  finite_index = find(~isinf(f));
  inf_index = find(isinf(f));
  xy = x0(nonz_indx,:);
  r = sqrt(sum(xy.^2,1));
  rnz = phi > 0;
  r_hat(:,rnz) = xy(:,rnz) ./ repmat(r(rnz),d-1,1);
  r_hat(:,~rnz) = repmat([0; 1],1,sum(~rnz)); % "r_hat" for axial point
  if (d == 3)
    e_theta = [r_hat(2,:); -r_hat(1,:)];
  end
  for i = 1:n_phi
    n_rays = size(dxinf{i},2);
    % Calculate projections of infinite-conjugate rays  (onto rhat, zhat, etheta)
    tanv_proj_inf(1,:) = sum(dxinf{i}(nonz_indx,:) .* repmat(r_hat(:,i),1,n_rays),1);
    zhatproj = dxinf{i}(z_indx,:);
    if (d == 3)
      tanv_proj_inf(2,:) = sum(dxinf{i}(nonz_indx,:) .* repmat(e_theta(:,i),1,n_rays),1);
    end
    % Calculate projections of finite-conjugate rays (onto t0, etheta)
    r_hat_d(nonz_indx,:) = r_hat(:,i);
    r_hat_d(z_indx,:) = 0;
    z_hat_d(z_indx,:) = 1;
    alpha_fc = lens.alpha{finite_index}(phi(i));
    t0 = r_hat_d .* repmat(sin(alpha_fc),d,1) + ...
      z_hat_d .* repmat(cos(alpha_fc),d,1);
    sN = r_hat_d * sin(phi(i)) + z_hat_d * cos(phi(i));
    ds = repmat(sN,1,n_rays) - sfin{i};
    tanv_proj_finite(1,:) = sum(ds .* repmat(t0,1,n_rays),1);
    if (d == 3)
      tanv_proj_finite(2,:) = sum(ds .* repmat([e_theta(:,i);0],1,n_rays),1);
    end
    % Calculate the "known" part of the coefficient relating them
    n = [lens.nL lens.nR];
    coef = n(finite_index) * lens.F{finite_index}(phi(i)) / n(inf_index);
    % Do the optimization
    Gguess = sum(coef * tanv_proj_finite(2,:) .* tanv_proj_inf(2,:)) ./ ...
	     sum(tanv_proj_inf(2,:).^2);
    p = [Gguess,pi/2];
    errfun = @(p) inf_err(p,coef*tanv_proj_finite,tanv_proj_inf,zhatproj);
    p = lsqnonlin(errfun,p);
    G(i) = p(1);
    alpha(i) = p(2);
  end
  % Load up the output lens structure
  lens_out = lens;
  
  lens_out.alpha{inf_index} = ppcreate(phi,alpha,'spline');
  lens_out.G{inf_index} = ppcreate(phi,G,'spline');
  
function err = inf_err(p,rhs,tanv_proj_inf,zhatproj)
  G = p(1);
  alpha = p(2);
  % change the sign in front of the cos to - so we get the tangent
  % direction for above-the-axis curvature
  err = (G * [sin(alpha)*tanv_proj_inf(1,:) - cos(alpha)*zhatproj; ...
    tanv_proj_inf(2,:)]) - rhs;
  fprintf('p: %g %g, total error: %g\n',p(1),p(2),sum(err(:).^2));
  