function phi = calcphi_single_path_series(image_series,H0,phi0,v)
% image_series: h-by-w-by-K
% v: actuator voltage for each image in series (a vector of length K)
% H0: pupil mask
% phi0: h-by-w-by-2, the first is the "constant" component and the second
%   is the "voltage-dependent" component,
%        phik = phic + v_k*phiv
%   This is the initial guess. The return will be the optimized version.
  
  %% Argument parsing
  impsz = size(image_series);
  if (length(impsz) ~= 3)
    error('Image data must be a 3-dimensional array')
  end
  imsz = impsz(1:2);
  K = impsz(3);
  
  % Initial guess
  if (nargin < 3 || isempty(phi0))
    % No guess supplied, make it flat. This can be a risky choice,
    % however, because of symmetries, so it's better to supply a guess.
    phi0 = zeros([imsz 2]);
  end

  % Create the function to be optimized
  optfunc = @(p) csps_pdwrapper(p,image_series,H0,v);
  
  % Do the optimization
  ops = struct('Display',true);
  ops.mu_max = 1;
  phi = conjgrad(optfunc,phi0,ops);
end

function [val,grad] = csps_pdwrapper(p,image_series,H0,v)
  impsz = size(image_series);
  K = impsz(3);
  phik = zeros(impsz);
  for k = 1:K
    phik(:,:,k) = p(:,:,1) + v(k)*p(:,:,2);
  end
  [val,gradphi] = pdpenalty(phik,image_series,H0);
  % Use the chain rule to get the derivatives with respect to phic and phiv
  grad = zeros([impsz(1:2) 2]);
  grad(:,:,1) = sum(gradphi,3);
  for k = 1:K
    grad(:,:,2) = grad(:,:,2) + v(k) * gradphi(:,:,k);
  end
  % Note! The next two lines mean that the derivative isn't right, but it
  % will make mu=1 be a reasonable size step
  mxgrad = max(abs(grad(:)));
  grad = grad / mxgrad;
end
