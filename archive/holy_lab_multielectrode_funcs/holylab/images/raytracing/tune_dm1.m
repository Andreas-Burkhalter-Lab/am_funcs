function DMp = tune_dm1(params,rays,virtual_source)
  infcorrect = false;
  if (nargin < 3)
    infcorrect = true;
  end
  
  if ~iscell(rays)
    rays = {rays};
  end
  n_groups = length(rays);
  
  thresh_rel = 1e-6;
  thresh_abs = 1e-6;
  n_iter_max = 100;

  n_vertices = length(params.DMparams);
  parallel = [cos(params.theta);sin(params.theta)];
  P = diag(parallel);
  normal = [parallel(2); -parallel(1)];
  params.normal = normal;
  lgrid = linspace(-params.length/2,params.length/2,n_vertices);
  
  if ~infcorrect
    if (virtual_source.normal(1) ~= 0)
      vsang = atan(virtual_source.normal(2)/virtual_source.normal(1));
    else
      vsang = pi/2;
      if (virtual_source.normal(2) < 0)
        vsang = -pi/2;
      end
    end
    vsN = diag(virtual_source.normal);
    vs_parallel = [-virtual_source.normal(2); virtual_source.normal(1)];
    vsP = diag(vs_parallel);
  end

  % Preprocess the rays: get their strike positions & angles
  % We can also compute the Hessian, which we approximate as constant
  % even for the finite-conjugate correction case using a paraxial
  % approximation (for the Hessian only, not for the gradient)
  A = zeros(n_vertices,n_vertices);
  if ~infcorrect
    zstrike = cell(1,n_groups);
    ystrike = cell(1,n_groups);
    ang_tau = cell(1,n_groups);
    dtheta_coef = cell(1,n_groups);
  else
    b = zeros(n_vertices,1);
  end
  for groupIndex = 1:n_groups
    % Calculate strike positions of the rays
    r = rays{groupIndex};
    r = r([r.valid]);
    n_rays = length(r);
    e = [r.e];
    ang = atan(e(2,:)./e(1,:));
    [t,strike_pos] = ray_flat_intersection(params,r);
    strike_dx = strike_pos - repmat(params.center,1,n_rays);
    % Convert this to a strike index (the mirror segment hit) and the
    % fractional offset (for linear interpolation)
    l = sum(P*strike_dx);
    segment = (l/params.length + 0.5)*(n_vertices-1);
    segmentIndex = floor(segment);
    segmentFrac = segment-segmentIndex;
    segmentIndex = segmentIndex+1; % unit offset
    keepFlag = (segmentIndex > 0 & segmentIndex < n_vertices);
    segmentIndex = segmentIndex(keepFlag);
    segmentFrac = segmentFrac(keepFlag);
    ang = ang(keepFlag);
    n_rays = sum(keepFlag);
    if infcorrect
      % In this case it's a true linear equation, so we can calculate
      % both sides of the equation
      a = sparse([1:n_rays 1:n_rays],[segmentIndex segmentIndex+1],...
        [1-segmentFrac segmentFrac],n_rays,n_vertices);
      bb = ang/2;
      bb = bb(:) - mean(bb);  % This is the "freedom to choose target"
    else
      % Here it's a nonlinear minimization, but we can approximate the
      % hessian as constant
      z = sum(vsN*strike_pos) - virtual_source.c;
      z = z(keepFlag);
      zstrike{groupIndex} = z;
      ytmp = sum(vsP*strike_pos);
      ystrike{groupIndex} = ytmp(keepFlag);
      a = sparse([1:n_rays 1:n_rays],[segmentIndex segmentIndex+1],...
        [z.*segmentFrac z.*(1-segmentFrac)],n_rays,n_vertices);
      % Skip the gradient calcuation for now, but do precompuations on
      % ray angles
      dtheta_coeff{groupIndex} = sparse([segmentIndex segmentIndex+1],...
        [1:n_rays 1:n_rays],...
        [1-segmentFrac segmentFrac],...
        n_vertices,n_rays);
      ang_tau{groupIndex} = ang - 2*params.theta + vsang;
    end
    a = full(a);
    a = a - repmat(mean(a,1),n_rays,1);
    A = A + a'*a;
    if infcorrect
      b = b + a'*bb;
    end
  end
  % Prepare to solve the linear equation. Since there may not be a ray
  % for each segment, the hessian may be singular. Therefore, we have to
  % use something like SVD.
  [U,S,V] = svd(A);
  s = diag(S);
  smax = max(abs(s));
  goodFlag = abs(s) > smax/1e6;
  sinv = zeros(1,length(s));
  sinv(goodFlag) = 1./s(goodFlag);
  Ainv = V*diag(sinv)*U';
  if infcorrect
    % We can just go ahead and solve it
    DMp = Ainv*b;
    return
  end
  
  % Now all that is left is the finite-conjugate case, which is
  % nonlinear. Do iterative optimization.
  is_done = false;
  DMp = zeros(n_vertices,1);
  iter = 0;
  while ~is_done
    grad = zeros(size(DMp));
    % Calculate the error & gradient
    err = 0;
    for groupIndex = 1:n_groups
      dtheta = 2*(DMp'*dtheta_coeff{groupIndex}) - ang_tau{groupIndex};
      n_rays = length(dtheta);
      yeff = zstrike{groupIndex} .* tan(dtheta);
      dy = ystrike{groupIndex} - mean(ystrike{groupIndex});
      ray_err = yeff - mean(yeff) - dy; % this is the current failure-to-focus
      err = err + sum(ray_err.^2);
      dre = repmat(2*(zstrike{groupIndex} ./ cos(dtheta).^2),n_vertices,1) ...
        .* dtheta_coeff{groupIndex};
      drec = dre - repmat(mean(dre,2),1,n_rays);
      grad = grad + 2*sum(repmat(ray_err,n_vertices,1) .* drec,2);
    end
    fprintf('%g...',err);
    % Calculate the correction
    deltaDMp = -Ainv * grad;
    % See whether the correction is an improvement
    err_old = err;
    lambda = 1;
    is_better = false;
    while (~is_better && lambda > 1e-6)
      DMptemp = DMp + lambda*deltaDMp;
      err = 0;
      for groupIndex = 1:n_groups
        dtheta = 2*(DMptemp'*dtheta_coeff{groupIndex}) - ang_tau{groupIndex};
        yeff = zstrike{groupIndex} .* tan(dtheta);
        dy = ystrike{groupIndex} - mean(ystrike{groupIndex});
        ray_err = yeff - mean(yeff) - dy;
        err = err + sum(ray_err.^2);
      end
      if (err > err_old)
        lambda = lambda/3;
      else
        is_better = true;
      end
    end
    if ~is_better
      error('Something is wrong')
    end
    DMp = DMptemp;
    iter = iter+1;
    is_done = (err_old - err < thresh_rel*(err+err_old) || err < thresh_abs || iter > n_iter_max);
  end
  fprintf('\n');