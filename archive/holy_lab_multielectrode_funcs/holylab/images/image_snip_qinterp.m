function [keepflagc,snip,snipg] = image_snip_qinterp(szout,offset,I,have_Icoef)
% image_snip_qinterp: cut out rectangular regions of an image, with subpixel interpolation
%
% Syntaxes:
%   [keepflagc,Isnip] = image_snip_qinterp(szout,offset,I)
%   [keepflagc,Isnip,Isnipg] = image_snip_qinterp(szout,offset,I)
% where
%    szout is a vector of length n_dims, giving the dimensions of the
%      output image (Isnip)
%    offset is a vector of length n_dims, specifying the coordinates of the
%      upper left corner of the region. If you have coordinates of the
%      center, then
%         offset = center - (szout-1)/2
%    I is the raw image. Alternatively you can supply the quadratic
%      interpolation coefficients, see the syntax below.
% On output,
%    keepflagc is a cell array of length n_dims, each element a logical
%      vector specifying whether the given "column" along that coordinate
%      is valid (if false, indicates the coordinate was beyond the edge)
%    Isnip is the snipped-out image. It may be smaller than requested, if
%      some of the coordinates were over the edge. If necessary, you can
%      pad it to full size in the following way:
%         Ipad = nan(szout);
%         Ipad(keepflagc{:}) = Isnip;
%    Isnipg (optional) is a cell array, with each element holding the
%      gradient of the image with respect to the corresponding coordinate.
%
% If you are snipping out many regions successively, it is to your
% advantage to pre-compute the interpolation coefficients like this:
%    Icoef = qinterp_grid_inverse(I,'nan');
% and then use the following syntax:
%   [keepflagc,...] = image_snip_qinterp(szout,offset,Icoef,true)
%
% See also: qinterp_grid_inverse.

% Copyright 2011 by Timothy E. Holy
  
  calc_grad = nargout > 2;
  if (nargin < 4)
    have_Icoef = false;
  end
  n_dims = length(szout);
  if iscell(szout)
    src = szout;
    szout = cellfun(@length,src);
    offset = offset + cellfun(@(x) x(1),src) - 1;
  end
  src = cell(1,n_dims);
  keepflagc = cell(1,n_dims);
  
  if all(offset == round(offset)) && ~calc_grad && ~have_Icoeff
    % Integer offset, just snip it out
    for dimIndex = 1:n_dims
      thissrc = (1:szout(dimIndex))+offset(dimIndex);
      keepflag = thissrc > 0 & thissrc <= size(I,dimIndex);
      src{dimIndex} = thissrc(keepflag);
      keepflagc{dimIndex} = keepflag;
    end
    snip = I(src{:});
  else
    % Sub-pixel offset, use quadratic interpolation
    if calc_grad
      [z,coef,coefg] = qinterp_coef(offset);
    else
      [z,coef] = qinterp_coef(offset);
    end
    offset = round(offset); % now just use the integer portion
    for dimIndex = 1:n_dims
      thissrc = (1:szout(dimIndex))+offset(dimIndex);
      % Need a buffer of 1 on all sides
      keepflag = thissrc > 1 & thissrc < size(I,dimIndex);
      src{dimIndex} = thissrc(keepflag);
      keepflagc{dimIndex} = keepflag;
    end
    snip = zeros(cellfun(@length,src),class(I));
    if calc_grad
      snipg = cell(1,n_dims);
      for dimIndex = 1:n_dims
        snipg{dimIndex} = sniptmp;
      end
    end
    if have_Icoef
      Icoef = I;
    else
      Icoef = qinterp_grid_inverse(I,'nan');
    end
    for zindx = 1:size(z,2)
      srctmp = src;
      for dimIndex = 1:n_dims
        srctmp{dimIndex} = srctmp{dimIndex}+z(dimIndex,zindx);
      end
      snip = snip + coef(zindx)*Icoef(srctmp{:});
      if calc_grad
        for dimIndex = 1:n_dims
          snipg{dimIndex} = snipg{dimIndex} + coefg(zindx)*Icoef(srctmp{:});
        end
      end
    end
  end
  