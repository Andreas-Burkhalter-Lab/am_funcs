function [p,ymean,yconflict] = flow_bn_stats(yfinal,l2final,w)
% FLOW_BN_STATS: statistics on convergence locations under balanced-neighborhood mean shift
% Syntax:
%   [p,ymean,ystd,lmean] = flow_bn_stats(yfinal,l2final,w)
% where
%   yfinal, l2final, and w are the outputs of FLOW_BN_WEIGHTS
% and
%   p is the probability assigned to each template/point combination
%   ymean is the mean convergence location of each point
%   yconflict is essentially a Z-score assigned to each data point,
%     indicating the uncertainty about convergence location on the scale of
%     the mean peak widths (i.e., the mean l2).
%
% See also: FLOW_BN_WEIGHTS.

% Copyright 2009 by Timothy E. Holy

  K = size(w,1);
  wsum = full(sum(w,1));
  p = w;
  for k = 1:K
    if issparse(w)
      j = find(p(k,:));
      ptmp = w(k,j) ./ wsum(j);
      p(k,j) = ptmp;
    else
      ptmp = w(k,:) ./ wsum;
      p(k,:) = ptmp;
    end
  end
  ymean = yfinal * p;
  yvar = yfinal.^2 * p;
  yvar = yvar - ymean.^2;
  l2mean = l2final * p;
  yconflict = sqrt(sum(yvar ./ l2mean,1));
end