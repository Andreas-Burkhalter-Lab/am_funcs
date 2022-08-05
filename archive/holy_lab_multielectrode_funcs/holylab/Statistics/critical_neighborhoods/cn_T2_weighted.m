function [val,g] = cn_T2_weighted(dx,w)
% cn_T2_weighted: compute the value of T2, and its gradient, as a function of weights
% dx = x - b, where b is the basepoint

  W = sum(w);
  X = sum(bsxfun(@times,w,dx),2);
  dX = bsxfun(@minus,X,W*dx);
  dX2 = sum(dX.^2,1);
  denom = sum(w.*dX2);
  val = W^2 * sum(X.^2) / denom;
  if (nargout > 1)
    dp = sum(bsxfun(@times,X,dx),1);  % dot product of X and dx
    wdX = sum(bsxfun(@times,dX,w),2);
    gdenom = dX2 + 2*sum(bsxfun(@times,wdX,dx),1) - 2*sum(X.^2) + 2*W*sum(w.*sum(dx.^2,1));
    g = 2*val/W  + (2*W^2/denom)*dp - (val/denom) * gdenom;
  end
end
