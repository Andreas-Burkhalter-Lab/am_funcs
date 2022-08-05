function [DMp,err] = tune_dm1_linear(params,rays,virtual_source)
  infcorrect = false;
  if (nargin < 3)
    infcorrect = true;
  end
  
  if ~iscell(rays)
    rays = {rays};
  end
  n_groups = length(rays);
  
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

  A = zeros(n_vertices,n_vertices);
  b = zeros(n_vertices,1);
  c = 0;
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
    % Set up the linear equation
    if infcorrect
      a = sparse([1:n_rays 1:n_rays],[segmentIndex segmentIndex+1],...
        [1-segmentFrac segmentFrac],n_rays,n_vertices);
      bb = ang/2;
    else
      z = sum(vsN*strike_pos) - virtual_source.c;
      a = sparse([1:n_rays 1:n_rays],[segmentIndex segmentIndex+1],...
        [z.*segmentFrac z.*(1-segmentFrac)],n_rays,n_vertices);
      y = sum(vsP*strike_pos);
      y0 = mean(y-z.*tan(2*params.theta - vsang - ang));
      unc_err = z.*tan(2*params.theta - vsang - ang) - (y-y0);
      %bb = y/2 + z.*ang/2;
      bb = y/2 - z.*tan(2*params.theta - vsang - ang)/2;
    end
    bb = bb(:) - mean(bb);
    c = c + sum(bb.^2);
    %bb = bb(:);
    a = full(a);
    a = a - repmat(mean(a,1),n_rays,1);
    A = A + a'*a;
    b = b + a'*bb;
  end
  % Since A might be singular, but it's square, we have to do the division
  % using SVD
  %DMp = A\b;
  [U,S,V] = svd(A);
  s = diag(S);
  smax = max(abs(s));
  goodFlag = abs(s) > smax/1e6;
  sinv = zeros(1,length(s));
  sinv(goodFlag) = 1./s(goodFlag);
  Ainv = V*diag(sinv)*U';
  DMp = Ainv*b;
%  DMp = DMp(end:-1:1);
  %DMp(segmentIndex)'.*(1-segmentFrac) + DMp(segmentIndex+1)'.*segmentFrac
  %a*DMp
  if (nargout > 1)
    err = DMp'*A*DMp - 2*(b'*DMp) + c;
  end
  
