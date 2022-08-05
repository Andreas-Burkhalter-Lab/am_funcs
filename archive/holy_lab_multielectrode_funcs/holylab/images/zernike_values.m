function Zval = zernike_values(rho,theta,coefList,normalizeStr)
% ZERNIKE_VALUES: compute the Zernike polynomials within a pupil
% Syntax:
%   Zval = zernike_values(rho,theta,coefList,normalizeStr)
% where
%   rho is the radius (pixels with rho > 1 or nan are assumed to lie
%     outside the pupil, and will have Zvalues of 0). This can be supplied
%     as a 2-dimensional array, so that the Zernike values are already
%     shaped to the pupil.
%   theta is the angle (its size and shape must match that of rho)
%   coefList is a list of Zernike polynomial indices. If supplied as a
%     1-by-K vector, we use the labeling scheme of zernfun2; if supplied
%     as a 2-by-K matrix, we use the labeling scheme of zernfun, where
%     the top row is n (the major number) and the bottom row is m (the
%     minor number).
%   normalizeStr (default 'norm'): has 3 acceptable values,
%     '': Zernike functions will be unnormalized
%     'norm': Zernike functions will be normalized in the continuous
%       sense, i.e., their values will be the analytical formulas
%       for normalized Zernike polynomials on grid locations
%     'normdiscrete': Zernike functions will be normalized on the grid,
%       so that sum(Z(:).^2) = N, the total number of pixels inside the
%       pupil.  As N tends to infinity, the results from 'normdiscrete'
%       tend to those from 'norm'.
% and
%   Zval is a [size(rho) K] array of Zernike polynomials with the desired
%     normalization scheme.

% Copyright 2008 by Timothy E. Holy
% Note: normalization scheme changed June 1 2009.

  %% Input parsing
  if (nargin < 4)
    normalizeStr = 'norm';
  end
  normalizeMode = strmatch(lower(normalizeStr),{'','norm','normdiscrete'},'exact');
  if isempty(normalizeMode)
    error(['Normalization string ' normalizeStr ' not recognized']);
  end
  [numberingMode,n_coefs] = size(coefList);
  if (numberingMode > 2)
    error('The coefficient list should be a row vector or a 2-by-K matrix');
  end
  numberingMode = 3-numberingMode; % 1 = (n,m); 2 = single-index mode
  if ~isequal(size(rho),size(theta))
    error('rho and theta are not of the same size');
  end
  
  %% Pupil boundaries
  rhosz = size(rho);
  H = zeros(rhosz,class(rho));
  pupilFlag = rho <= 1;
  H(pupilFlag) = 1;
  normPupil = sum(H(:));
  if (normPupil == 0)
    error('There are no pixels in your pupil!');
  end
  
  %% Calculate the values
  args = {rho(pupilFlag),theta(pupilFlag)};
  if (normalizeMode == 2)
    args{3} = 'norm';
  end
  if (numberingMode == 1)
    Ztmp = zernfun(coefList(1,:),coefList(2,:),args{:});
  else
    Ztmp = zernfun2(coefList,args{:});
  end
  Zval = zeros([rhosz n_coefs],class(rho));
  for indx = 1:n_coefs
    tmp = H;
    tmp(pupilFlag) = Ztmp(:,indx);
    if (normalizeMode == 3)
      tmp = tmp / sqrt(sum(Ztmp(:,indx).^2)/normPupil); % normalize on discrete grid
    end
    Zval(:,:,indx) = tmp;
  end
end