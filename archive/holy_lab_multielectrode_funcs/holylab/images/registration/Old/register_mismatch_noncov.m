function [imMg,val,grad] = register_mismatch_noncov(g,imM,imF,theta)
% REGISTER_MISMATCH_NONCOV: compute the mismatch error under deformation
%
% This function implements warping and penalty-function evaluation for
% non-covariant image registration. The penalty is
%      E = sum(theta(g(x)) .* dI(x).^2)
% where dI = imFixed - imMoving(g(x)) and theta is a weight function.
% Usually, theta should be zero along the edge of the image (or neighboring
% NaNs) and 1 everywhere else.
%
% Syntax:
%   imMg = register_mismatch_noncov(g,imM)
% This version simply evaluates the moving image at the positions defined
% by the array g (with size [size(imMg) ndims(imM)]).
%
%   [imMg,val] = register_mismatch_noncov(g,imM,imF,theta)
% This version also evaluates the penalty, supplying the fixed image imF
% and the weight function theta.
%
%   [imMg,val,grad] = register_mismatch_noncov(g,imM,imF,theta)
% This version computes the gradient of the penalty with respect to g.

% Copyright 2009-2010 by Timothy E. Holy
  
  normalize = true;  % controls whether we divide by the # of _valid_ pixels

  szM = size(imM);
  dimFlag = szM > 1;
  if ~dimFlag(1)
    error('First dimension must not be unity (try supplying the transpose image)');
  end
  n_dims = length(szM(dimFlag));
  
  if (nargout > 1)
    wF = ~isnan(imF);
  end
  
  %% This version uses the MEX file imqinterp---comment out if you don't
  % want to use it
  if (nargout < 2)
    imMg = imqinterp(g,imM);
    return
  else
    if (nargout < 3)
      out = cell(1,2);
    else
      out = cell(1,4);
    end
    [out{:}] = imqinterp(g,imM,theta);
    imMg = out{1};
    wMg = out{2};
    dI = imF - imMg;
    wMg(isnan(wMg)) = 0;
    dI(isnan(dI)) = 0;
    wprod = wMg.*wF;
    if normalize
      W = sum(wprod(:));
    else
      W = numel(imF);
    end
    val = sum(wprod(:) .* dI(:).^2) / W;
    if (nargout > 2)
      dimMdg = out{3};
      dwMdg = out{4};
      repsz = [ones(1,n_dims) n_dims];
      gradm = dwMdg.*repmat(wF.*dI.^2/W,repsz)...
        - repmat((2/W)*wprod.*dI,repsz).*dimMdg;
      if normalize
        gradm = gradm - repmat((val/W)*wF,repsz).*dwMdg;
      end
      gradm(isnan(gradm)) = 0;
%       szc = num2cell(szM(dimFlag),1);
%       grad = mat2cell(gradm,szc{:},ones(1,n_dims));
      grad = gradm;
    end
  end
  return
      

  if ~iscell(gc)
    error('This function expects cell-array g');
  end
  gc = gc(szM > 1);
  szM1 = szM(szM > 1);
  n_dims = length(gc);
  
  %% Prepare g and get all the indexing done
  gc_round = cell(size(gc));
  dg = cell(size(gc));
  for i = 1:n_dims
    gc_round{i} = round(gc{i});
    dg{i} = gc{i} - gc_round{i};
    gc_round{i}(gc_round{i} < 1 | gc_round{i} > szM1(i)) = nan;
  end
  if (n_dims > 1)
    ind = sub2ind(szM1,gc_round{:});
  else
    ind = gc_round{1};
  end
  dimoffset = [1 cumprod(szM1)];
  z = zeros(3^n_dims,n_dims);
  memoffset = zeros(1,3^n_dims);
  % z will be a counter running from -1 to 1, and we have to increment
  % with carrying
  z(1,:) = -1;
  for i = 1:3^n_dims
    memoffset(i) = sum(z(i,:) .* dimoffset(1:end-1));
    znext = z(i,:);
    znext(1) = znext(1) + 1;
    j = 1;
    while (znext(j) > 1)
      % the carry operation
      znext(j) = -1;
      if (j < n_dims)
        j = j+1;
        znext(j) = znext(j)+1;
      end
    end
    if (i < size(z,1))
      z(i+1,:) = znext;
    end
  end
  
  %% Create the warped image and warped weight
  if (nargout < 2)
    imMg = warp_im(ind,memoffset,z,dg,gc_round,imM);
    return
  end
  [imMg,w] = warp_im(ind,memoffset,z,dg,gc_round,imM,theta);
  w(isnan(w)) = 0;
  if any(isnan(imMg(:)) & w(:) ~= 0)
    error('Something weird happened');
  end
  
  %% Calculate the mismatch penalty
  dI = imF - imMg;
  if normalize
    W = sum(w(:));
  else
    W = numel(imF);
  end
  val = nansum(w(:) .* dI(:).^2) / W;
  if (nargout < 3)
    return
  end
  
  %% Calculate the terms in the gradient
  [dI2dg,dwdg] = grad_warp_im(ind,memoffset,z,dg,gc_round,imM,theta);
  grad = cell(1,n_dims);
  for i = 1:n_dims
    if normalize
      grad{i} = dwdg{i}/W.*(dI.^2 - val) - (2/W)*w.*dI.*dI2dg{i};
    else
      grad{i} = dwdg{i}/W.*dI.^2 - (2/W)*w.*dI.*dI2dg{i};
    end
    grad{i}(isnan(grad{i})) = 0;
  end
  if (length(szM) > length(szM1))
    gradtmp = grad;
    grad = cell(1,length(szM));
    grad(szM > 1) = gradtmp;
    for i = 1:length(szM)
      if isempty(grad{i})
        grad{i} = zeros(szM,class(imM));
      end
    end
  end

  
function varargout = warp_im(ind,memoffset,z,dg,gc_round,varargin)
  n_dims = length(dg);
  szg = size(dg{1});
  goodFlag = cell(1,length(varargin));  % indicates whether pixels need over-the-edge values
  imtmp = cell(1,length(varargin)); % holds results for valid pixels in a column vector
  for imIndex = 1:length(varargin)
    goodFlag{imIndex} = true(szg);
    for dimIndex = 1:n_dims
      % Check to see if gc_round will require an over-the-edge pixel
      goodFlag{imIndex} = goodFlag{imIndex} & gc_round{dimIndex} > 1 & gc_round{dimIndex} < size(varargin{imIndex},dimIndex);
    end
    varargout{imIndex} = zeros(szg,class(varargin{imIndex}));
    varargout{imIndex}(~goodFlag{imIndex}) = nan;
    imtmp{imIndex} = zeros(sum(goodFlag{imIndex}(:)),1,class(varargin{imIndex}));
  end
  for i = 1:size(z,1)
    coef = ones(szg,class(dg{1}));
    for j = 1:n_dims
      if z(i,j)
        coef = coef .* (dg{j}+ z(i,j)/2).^2/2;
      else
        coef = coef .* (0.75-dg{j}.^2);
      end
    end
    for imIndex = 1:length(varargin)
      cind = ind(goodFlag{imIndex}) + memoffset(i);
      imtmp{imIndex} = imtmp{imIndex} + coef(goodFlag{imIndex}).*varargin{imIndex}(cind);
    end
  end
  for imIndex = 1:length(varargin)
    varargout{imIndex}(goodFlag{imIndex}) = imtmp{imIndex};
  end

  
function varargout = grad_warp_im(ind,memoffset,z,dg,gc_round,varargin)
  n_dims = length(dg);
  szg = size(dg{1});
  goodFlag = cell(1,length(varargin));
  for imIndex = 1:length(varargin)
    goodFlag{imIndex} = true(szg);
    varargout{imIndex} = cell(1,n_dims);
    for dimIndex = 1:n_dims
      varargout{imIndex}{dimIndex} = zeros(szg,class(varargin{imIndex}));
    end
  end
  for dimIndex = 1:n_dims
    for i = 1:size(z,1)
      coef = ones(szg,class(dg{1}));
      for j = 1:n_dims
        if (j == dimIndex)
          if z(i,j)
            coef = coef .* (dg{j}+z(i,j)/2);
          else
            coef = coef .* (-2*dg{j});
          end
        else
          if z(i,j)
            coef = coef .* (dg{j}+ z(i,j)/2).^2/2;
          else
            coef = coef .* (0.75-dg{j}.^2);
          end
        end
        % Check each pixel to see if it's within the domain of each image
        for imIndex = 1:length(varargin)
          goodFlag{imIndex} = goodFlag{imIndex} & gc_round{j} > -z(i,j) & gc_round{j} <= size(varargin{imIndex},j)-z(i,j);
        end
      end
      cind = ind + memoffset(i);
      for imIndex = 1:length(varargin)
        varargout{imIndex}{dimIndex}(goodFlag{imIndex}) = varargout{imIndex}{dimIndex}(goodFlag{imIndex}) + coef(goodFlag{imIndex}).*varargin{imIndex}(cind(goodFlag{imIndex}));
      end
    end
    for imIndex = 1:length(varargin)
      varargout{imIndex}{dimIndex}(~goodFlag{imIndex}) = nan;
    end
  end

  