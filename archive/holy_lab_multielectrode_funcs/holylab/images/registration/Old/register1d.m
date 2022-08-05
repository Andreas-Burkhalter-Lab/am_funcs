function [g,img,Aout,x0out] = register1d(im1,im2,g0,lambda,A,x0,options)
  x2 = 1:length(im2);
  if (lambda < inf)
    % g0 is the input g(x1), where x1 are the coords of im1
    if isempty(g0)
      g0 = 1:length(im1);  % the identity transformation
    end
    img = register1d_newimg(x2,im2,g0);
    lKeep = ~isnan(img);  % Contains the pixels within the support of both
    gradim1 = gradient(im1(lKeep));
    gradimg = gradient(img(lKeep));
    w = img(lKeep).*gradim1 - im1(lKeep).*gradimg;
    wcs = cumsum(w);  % the integral of w
    J_base = exp(wcs/lambda);
    % Extrapolate J into regions where it is not defined
    % (want to use nearest neighbor, basically to continue J)
    J_base = interp1(find(lKeep),J_base,1:length(im1),'nearest','extrap');
  else
    J_base = ones(1,length(im1));
    lKeep = true(size(J_base));
  end
  g_base = cumsum(J_base);
  g = nan(size(lKeep));
  if ~options.optaffine
    g = A*g_base + x0;
  else
    p = fminsearch(@(p) register1d_mse(p,x2,im1,im2,g_base),[A x0]);
    g = p(1)*g_base + p(2);
  end
  if (nargout > 1)
    img = register1d_newimg(x2,im2,g);
    Aout = p(1);
    x0out = p(2);
  end
  
function img = register1d_newimg(x2,im2,g)
  im2gx = interp1(x2,im2,g);
  J = gradient(g);
  img = sqrt(abs(J)) .* im2gx;  % This is the warped im2

function mse = register1d_mse(p,x2,im1,im2,g_base)
  g = p(1)*g_base + p(2);
  img = register1d_newimg(x2,im2,g);
  lKeep = ~isnan(img);
  mse = mean((im1(lKeep) - img(lKeep)).^2);
  