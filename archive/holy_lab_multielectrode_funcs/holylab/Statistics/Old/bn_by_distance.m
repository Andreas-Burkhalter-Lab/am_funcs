function [nbrInfo,sq_summed_disp,summed_sq_disp,options] = bn_by_distance(x,y,ell,options)
% BN_BY_DISTANCE: get the balanced neighborhood, ranking by distance
% Syntax:
%   nbrInfo = bn_by_distance(x,y)
%   nbrInfo = bn_by_distance(x,y,ell)
%   [nbrInfo,sq_summed_disp,summed_sq_disp] = bn_by_distance(x,y,ell,options)

  [d,N] = size(x);
  if (nargin < 3)
    ell = ones(d,1);
  end
  if (nargin < 4)
    options = struct;
  end
  options = default(options,'mode','discrete');
  
  dx = x - repmat(y,1,N);
  dx = dx ./ repmat(ell,1,N);
  if strcmp(options.mode,'gaussian')
    [w,f1,f2,kappa] = bn_gaussian(dx);
    % Calculate the relevant fields in nbrInfo
    wsum = sum(w);
    ynew = y + f1/wsum;
    ynew = ynew .* ell;
    l2new = l2 .* (f2/wsum);
    l2new(l2new == 0) = l2(l2new == 0);
    nbrInfo = struct('ybase',y,'ell',ell,'neighborDist',sqrt(d2),'kappa',kappa,'n',wsum,'n_z',wsum,'nbrhoodMean',ynew,'nbrhoodMSDisp',l2new,'RMSDist',sqrt(sum(d2.*w)/wsum));
    return
  end
  R2 = sum(dx.^2,1);
  [sR2,sortOrder] = sort(R2);
  [dxmean,dx2mean,n,n_exceed,sq_summed_disp,summed_sq_disp,options] = bn_preordered(dx(:,sortOrder),options);
  ynew = y + dxmean .* ell;
  l2new = dx2mean .* ell.^2;
  nbrInfo = struct('ybase',y,...
    'ell',ell,...
    'neighborLabel',sortOrder,...
    'neighborDist',sqrt(sR2),...
    'n',n,...
    'n_z',n_exceed,...
    'nbrhoodMean',ynew,...
    'nbrhoodMSDisp',l2new,...
    'RMSDist',sqrt(mean(sR2(1:n))));
  %nbrInfo.kurtosis = kurtosis(nbrInfo.neighborDist(1:n));
end
