function hline = plot_covar_as_ellipse(varargin)
% plot_covar_as_ellipse: plot covariance matrix as ellipse
% Syntax:
%   hline = plot_covar_as_ellipse(mu,C)
%   hline = plot_covar_as_ellipse(hax,mu,C)
% where
%   mu is a 2-vector giving the center of the ellipse
%   C is the 2-by-2 covariance matrix
% and
%   hline is the line handle for the ellipse
  
% Copyright 2010 by Timothy E. Holy
  
  if ishandle(varargin{1})
    hax = varargin{1};
    mu = varargin{2};
    C = varargin{3};
  else
    hax = gca;
    mu = varargin{1};
    C = varargin{2};
  end
  
  d = length(mu);
  if (d ~= 2)
    error('works only for 2 dimensions');
  end
  [R,p] = chol(C);
  if (p == 0)
    theta = linspace(0,2*pi,101);
    dx = R'*[cos(theta);sin(theta)];
    hline = line(mu(1)+dx(1,:),mu(2)+dx(2,:),'parent',hax);
  else
    % C is not positive definite, let's plot a line along the nonzero
    % eigenvalue
    [V,D] = eig(C);
    l = sqrt(D(2,2));
    dx = V(:,2)*[1 -1];
    hline = line(mu(1)+dx(1,:)*l/2,mu(2)+dx(2,:)*l/2);
  end
  