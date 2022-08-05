function clust = climb_gaussian(X,alpha,mu,sigma2,options)
% CLIMB_GAUSSIAN: cluster by moving points up a mixture-of-gaussians density
%
% The movement up the density is given by the meanshift algorithm, see
% Comaniciu & Meer. Their equations are generalized to allow the alpha
% pre-factor.
%
% Syntax:
%   clust = climb_gaussian(X,alpha,mu,sigma2,options)
% where
%   X is a d-by-N matrix of data points;
%   alpha is a ngaussians-by-1 vector of amplitudes;
%   mu is a d-by-ngaussians matrix of center positions;
%   sigma2 is a 1-by-ngaussians vector of widths;
%   options is a structure which may have the following fields:
%     tol_dX2: the tolerance on the change in X values, in terms of the
%       _square_ euclidean distance between x_i(old) and x_i(new)
%       (default: sqrt(eps)*mean(sigma2))
%     ploteach: if true, plots x on each iteration (default false);
%     plotpause: if true, pauses the plot on each iteration (default
%       false);
%     plotclust: if true, plots each cluster with colored points at the end
%       of clustering (default false);
%     plotproject: if supplied, projects the data & centers down to two
%       dimensions before plotting (a 2-by-d matrix) (default: use first
%       two dimensions of the input);
% and
%   clust is the cluster number assigned to each data point
%
% See also: MEANSHIFT.

  if (nargin < 5)
    options = struct;
  end
  [d,N] = size(X);
  options = climb_gaussian_options(options,sigma2,d);
  
  if options.ploteach
    climb_gaussian_plot(X,mu,options);
  end
  % Iteratively move points
  Xold = X;
  Xnew = climb_gaussian_iter(Xold,alpha,mu,sigma2);
  while any(sum((Xnew-Xold).^2,1) > options.tol_dX2)
    if options.ploteach
      climb_gaussian_plot(Xnew,mu,options);
    end
    Xold = Xnew;
    Xnew = climb_gaussian_iter(Xold,alpha,mu,sigma2);
  end
  % Points have (basically) stopped. Assign each point to its (now)
  % closest gaussian center; points that end up at the same center are in
  % the same cluster.
  [md,indexCenter] = mindist(Xnew,mu);
  [uindex,tmp,clust] = unique(indexCenter);
  if options.plotclust
    Xp = options.plotproject*X;
    clabel = agglabel(clust);
    nclust = length(clabel);
    cla
    marker = '.';
    if (N < 500)
      marker = 'o';
    end
    for i = 1:nclust
      col = unique_color(i,nclust);
      line(Xp(1,clabel{i}),Xp(2,clabel{i}),...
        'MarkerFaceColor',col,...
        'MarkerEdgeColor',col,...
        'LineStyle','none',...
        'Marker',marker);
    end
  end
  

function Xnew = climb_gaussian_iter(X,alpha,mu,sigma2)
  [d,N] = size(X);
  [d,q] = size(mu);
  sd = sqrdist(X,mu);
  sigma2rep = repmat(sigma2,N,1);
  coef = alpha./(2*pi*sigma2).^(d/2);
  sdsig = sd ./ sigma2rep;
  p = repmat(coef,N,1).*exp(-sdsig/2); % p(i,j) refers to x_i & mu_j
  pnorm = p ./ repmat(sum(p,2),1,q);  % fraction of x_i explained by j
  Xnew = mu*pnorm';

  
function climb_gaussian_plot(X,mu,options)
  mup = options.plotproject*mu;
  Xp = options.plotproject*X;
  plot(Xp(1,:),Xp(2,:),'.')
  hold on
  scatter(mup(1,:),mup(2,:),'ro');
  hold off
  axis equal
  if options.plotpause
    pause
  else
    drawnow
  end


function options = climb_gaussian_options(options,sigma2,d)
  if ~isfield(options,'ploteach')
    options.ploteach = 0;
  end
  if ~isfield(options,'plotpause')
    options.plotpause = 0;
  end
  if ~isfield(options,'plotclust')
    options.plotclust = 0;
  end
  if ~isfield(options,'plotproject')
    options.plotproject = [eye(2), zeros(2,d-2)];
  end
  if ~isfield(options,'tol_dX2')
    options.tol_dX2 = sqrt(eps)*mean(sigma2);
  end
  