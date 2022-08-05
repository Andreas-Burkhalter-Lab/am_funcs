function T2 = cn_T2(nvec,mu,C,covarianceModel,T2thresh)
  haveThresh = nargin > 4;
  switch covarianceModel
    case 'isotropic'
      T2 = nvec.*sum(mu.^2,1)./C;  % FIXME this does not perfectly agree with cn_preordered
    case 'diagonal'
      T2 = nvec.*sum(mu.^2./C,1);
    case 'full'
      [d,N] = size(mu);
      T2 = zeros(1,N);
      for i = d+1:N
        T2(i) = nvec(i)*(mu(:,i)'*(C(:,:,i)\mu(:,i)));
        if (haveThresh && T2(i) > T2thresh(i))
          break
        end
      end
    otherwise
      error('Not yet implemented')
  end
end
  