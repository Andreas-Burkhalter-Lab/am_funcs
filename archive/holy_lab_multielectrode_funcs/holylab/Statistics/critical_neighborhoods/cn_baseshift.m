function [dp,info] = cn_baseshift(n_all,p_all,C_all,params)
% cn_baseshift: update model parameters via critical neighborhoods

  %% Move towards first significant centroid displacement
  T2 = cn_T2(n_all,p_all,C_all,params.covarianceModel,params.T2thresh);
  info.T2 = T2;
  n = find(T2 > params.T2thresh,1,'first');
  if isempty(n)
    % Nothing violated significance, let's quit
    info.n = length(T2);
    dp = [];
    return
  end
  deltap1 = p_all(:,n);
  info.n = n-1;
  [alpha1,alphamax,~,info.alphaindex] = cn_linesearch(p_all,C_all,deltap1,params);
  %% Move towards centroid that most constrained the previous search
  p_all = bsxfun(@minus,p_all,alpha1*deltap1); % correct the centroid displacements for the new basepoint
  cutoff = find(alphamax < 0,1,'first');  % avoid searches in opposite direction
  if isempty(cutoff)
    cutoff = length(alphamax);
  else
    cutoff = cutoff-1;
  end
  cutoff = min(cutoff,info.alphaindex);
  [~,indx] = min(alphamax(1:cutoff));  % index of most constraining nbrhood
  deltap2 = p_all(:,indx);
  info.interiorIndex = indx;
  % Make the new search direction conjugate to the first
  deltap2 = deltap2 - proj(deltap2,deltap1,C_all,info.alphaindex,params.covarianceModel); % subtract component parallel to deltap1
  alpha2 = cn_linesearch(p_all(:,1:info.alphaindex),C_all,deltap2,params);
  if (alpha2 > 1)
    alpha2 = 1;  % don't move any farther than the constraining centroid
  end
  dp = alpha1*deltap1 + alpha2*deltap2;
  info.alpha = [alpha1 alpha2];
end

function vp = proj(v,u,C_all,index,covarianceModel)
  switch covarianceModel
    case 'isotropic'
      vp = sum(v.*u)/sum(u.^2)*u;
    case 'diagonal'
      C = C_all(:,index);
      vp = sum(v.*u./C)/sum(u.^2./C)*u;
    case 'full'
      C = C_all(:,:,index);
      Cinvu = C\u;
      vp = sum(v'*Cinvu)/sum(u'*Cinvu)*u;
  end
end