function phi = calcphi_split_path(image_pair,H0,phi0,roughness_penalty)
% image_pair: h-by-w-by-2 (first of pair is "unaberrated")
% pupil_data:
  
  %% Argument parsing
  impsz = size(image_pair);
  if (length(impsz) ~= 3)
    error('Image data must be a 3-dimensional array')
  end
  if (impsz(3) ~= 2)
    error('Must be supplied as image pairs');
  end
  if (nargin < 4)
    roughness_penalty = 0;
  end
  
  % Initial guess
  if (nargin < 3 || isempty(phi0))
    % No guess supplied, make it flat. This can be a risky choice,
    % however, because of symmetries, so it's better to supply a guess.
    phi0 = zeros(impsz(1:2));
  end

  % Create the function to be optimized
  optfunc = @(p) csp_pdwrapper(p,image_pair,H0,roughness_penalty);
  
  % Do the optimization
  %optoptions = optimset('Display','iter','GradObj','on');
  %p = fminunc(optfunc,p0,optoptions);
  ops = struct('Display',true);
%   if (roughness_penalty == 0)
    ops.mu_max = 1;
%   end
  phi = conjgrad(optfunc,phi0,ops);
end

function [val,grad] = csp_pdwrapper(p,image_pair,H0,roughness_penalty)
  impsz = size(image_pair);
  phi = zeros(impsz);
  phi(:,:,2) = p;
  [val,gradphi] = pdpenalty(phi,image_pair,H0);
  grad = gradphi(:,:,2);
  if (roughness_penalty ~= 0)
    ps = fftshift(p);
    H0s = fftshift(H0);
    llapp = roughness_penalty * del2(H0s.*ps);
    pllapp = ps .* llapp;
    val = val - (1/2) * sum(pllapp(:));
    grad = grad - fftshift(llapp);
  end
%   else
    % Note! The next two lines mean that the derivative isn't right, but it
    % will make mu=1 be a reasonable size step
    mxgrad = max(abs(grad(:)));
    grad = grad / mxgrad;
%   end
end
