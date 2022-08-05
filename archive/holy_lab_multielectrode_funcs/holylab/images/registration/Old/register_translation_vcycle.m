function [dx,im2shift,err] = register_translation_vcycle(dx0,im2,options)
% REGISTER_TRANSLATION_VCYCLE: implement one v-cycle of shift registration
%
% See also: REGISTER_TRANSLATION_MULTIRES.

% Copyright 2007 by Timothy E. Holy

  output_progress = false;
  if isfield(options,'display')
    output_progress = options.display;
  end
  grids = options.image_grid;
  
  % Shift im2 to the best starting location
  if ~isfield(options,'im2shift')
    if all(dx0 == 0)
      im2tmp = im2;
    else
      im2tmp = image_shift(im2,dx0);
    end
  else
    im2tmp = options.im2shift;
  end
  % Create pyramid of lower-resolution versions
  grids(1).psiM = im2tmp;
  n_grids = length(grids);
  for gridIndex = 2:n_grids
    grids(gridIndex).psiM = imrestrict_mex(grids(gridIndex-1).psiM,...
				       grids(gridIndex).decimate);
  end
  % For the lowest-resolution image pair, find the shift that produces
  % the optimal alignment
  dx = 0*dx0;
  mu = 1;
  errOld = inf;
  %err = rtv_mismatch(grids(end).psiF,grids(end).psiM,dx,dx,0);
  err = imcompare(grids(end).psiF,{grids(end).psiM});
  if output_progress, fprintf('At coarsest grid\nInitial error = %g\n',err); end
  iter = 1;
  while ((errOld - err) > 1e-3*err || err > errOld)
    errOld = err;
    [dx,mu,err] = rtv_dec(grids(end).psiF,grids(end).psiM,dx,mu);
    if (iter == 1)
      mustart = mu;
    end
    if output_progress
      fprintf('iter %d: err = %g, mu = %g, dx = ',iter,err,mu);
      fprintf('%g ',dx);
      fprintf('\n');
    end
    mu = 2*mu;  % Default to a slightly larger mu next time
    iter = iter+1;
  end
  
  % Now work back up the hierarchy, applying the previous best-guess, 
  % but take only one descent step
  mu = mustart;
  for gridIndex = n_grids-1:-1:1
    scale_fac = grids(gridIndex+1).pixel_spacing ./ ...
      grids(gridIndex).pixel_spacing;
    dx = dx .* scale_fac;
    if output_progress
      fprintf('Scale %d: input dx = ',gridIndex)
      fprintf('%g ',dx);
    end
    if (gridIndex > 1)
      [dx,mu] = rtv_dec(grids(gridIndex).psiF,grids(gridIndex).psiM,dx,mu);
    else
      % For the highest-resolution, bias by the original offset and use
      % the entire input (that way we don't get more truncation due to
      % edge effects than absolutely necessary)
      [dx,mu,err,im2shift] = rtv_dec(grids(gridIndex).psiF,im2,dx+dx0,mu);
    end
    if output_progress
      fprintf('; output dx = ')
      fprintf('%g ',dx);
      fprintf(' (mu = %g)\n',mu);
    end
  end
  %[err,im2shift] = rtv_mismatch(grids(1).psiF,im2,dx,0,0);
  if output_progress
    fprintf('Cumulative dx: ');
    fprintf('%g ',dx);
    fprintf('\n');
  end

function [dx,mu,err,ims] = rtv_dec(im1,im2,dx0,mu)
  if (nargin < 4)
    mu = 1;
  end
  % Calculate the shifted image at the basepoint, and the gradient
  %R = makeresampler('linear','fill');
  %im2tmp = image_shift(im2,dx0,R);
  im2tmp = image_shift(im2,dx0);
  diffim = im1-im2tmp;
  sz = size(im2);
  dimFlag = sz > 1;
  dimIndex = find(dimFlag);
  n_dims = sum(dimFlag);
  step = zeros(1,n_dims);
  for dimIterator = 1:n_dims
    gradim = deriv(im2tmp,dimIndex(dimIterator));
    tmp = diffim .* gradim;
    step(dimIterator) = double(nanmean(tmp(:)));
  end
%  fmu = @(mu) rtv_mismatch(im1,im2,dx0,step,mu);
  fmu = @(mu) rtv_shiftray(im2,dx0,step,mu);
  %ops.mu = mu;
  %ops.bracketonly = true;
  %[mu,err] = lindec(fmu,ops);
  %[mu,err] = linmin(fmu,ops);
  ops = struct('f0',im2tmp,...
    'compare_func',@(imc) imcompare(im1,imc));
  [mu,ims,err] = lindip(fmu,mu,ops);
  dx = dx0 + step*mu;
  
function [err,im2tmp] = rtv_mismatch(im1,im2,base,step,mu)
  dx = base + double(mu)*step;
  if any(isinf(dx)) || any(isnan(dx))
    err = inf;
    return
  end
  %R = makeresampler('linear','fill');
  %im2tmp = image_shift(im2,dx,R);
  im2tmp = image_shift(im2,dx);
  diffim2 = (im1-im2tmp).^2;
  err = double(nanmean(diffim2(:)));
  if isnan(err)
    err = inf;
  end

function ims = rtv_shiftray(im,base,step,mu)
  dx = base + mu*step;
  if any(isinf(dx)) || any(isnan(dx))
    ims = nan(size(im),class(im));
    return
  end
  ims = image_shift(im,dx);
