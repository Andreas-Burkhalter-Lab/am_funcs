function [alpha,alphamax,alphamin,index] = cn_linesearch(mu,C,deltax,params)
  % mu is the neighborhood mean relative to the base point, e.g.,
  %      mean(x(:,nbrList),2) - x0
  % C is the neighborhood covariance
  % deltax is the search direction
  % params has fields
  %
  % alpha is the coefficient of deltax in the largest-possible step
  % alphamax is the list of upper constraints, in nbrhood order
  % alphamin is the list of lower constraints, in nbrhood order
  % index is the index of the last mutually-consistent set of constraints
  
  % Copyright 2011 by Timothy E. Holy
  
  [~,N] = size(mu);
  nvec = 1:N;
  
  %% Define the linesearch boundaries for each candidate neighborhood
  switch params.covarianceModel
    case 'isotropic'
      A = params.T2lsthresh(1:N).*C./nvec;
      a = sum(deltax.^2,1)./A;
      b = sum(bsxfun(@times,deltax,mu),1)./A;
      c = sum(mu.^2,1)./A - 1;
    case 'diagonal'
      A = bsxfun(@times,params.T2lsthresh(1:N)./nvec,C);
      a = sum(deltax.^2./A,1);
      b = sum(bsxfun(@times,deltax,mu./A),1);
      c = sum(mu.^2./A,1) - 1;
    case 'full'
      if size(C,3) > N
        C = C(:,:,1:N);
      end
      A = bsxfun(@times,reshape(params.T2lsthresh(1:N)./nvec,[1 1 N]),C);
      a = zeros(1,N);
      b = a; c = a;
      istart = find(~isinf(params.T2lsthresh(1:N)),1,'first');
      if isempty(istart)
        alpha = 1;
        alphamin = zeros(1,N);
        alphamax = inf(1,N);
        index = N;
        return
%         error('All values of T2 are inf');
      end
      for i = istart:N
        Aiterms = A(:,:,i)\[mu(:,i) deltax];
        a(i) = deltax'*Aiterms(:,2);
        b(i) = deltax'*Aiterms(:,1);
        c(i) = mu(:,i)'*Aiterms(:,1) - 1;
        if (~isinf(params.T2lsthresh(i)) && b(i)^2 - a(i)*c(i) < 0)
          break;  % ellipse miss, so we know we are done
        end
      end
    otherwise
      error('Not implemented')
  end
  % Find the first ellipse that is entirely missed by the search line
  % (we'll truncate consideration of neighborhoods at that point)
  q = b.^2 - a.*c;
  nmax = find(q < 0,1,'first');
  if isempty(nmax)
    nmax = N;
  else
    nmax = nmax-1;
  end
  % For each ellipse, calculate the range of alpha values that lie within the ellipse
  sqrtq = sqrt(q(1:nmax)); sqrtq(end+1:N) = 0;
%   b = b(1:nmax);
%   a = a(1:nmax);
  alphamin = (b-sqrtq)./a;
  alphamax = (b+sqrtq)./a;
  alphamax(isnan(alphamax)) = inf;
  if params.positive_only
    alphamin(isnan(alphamin)) = 0; % 0, rather than -inf, because we want to search alpha >= 0
  else
    alphamin(isnan(alphamin)) = -inf;
  end
  
  %% Find the largest neighborhood & step satisfying all constraints
  calphamin = cummax(alphamin);
  calphamax = cummin(alphamax);
  index = find(calphamin >= calphamax,1,'first');
  if isempty(index)
    index = nmax;
  else
    index = index-1;
  end
  if isinf(params.T2lsthresh(index)) || isnan(params.T2lsthresh(index))
    % There is no valid interval consistent with the current point, other
    % than those marked as non-restrictive due to nMin requirements
    % Here, our best bet is to go along deltax to the point closest to the
    % center of the first restrictive ellipse
    alpha = b(index+1)/a(index+1);
  else
%     alpha = (calphamax(index) + calphamin(index))/2;  % midpoint of last valid interval
    f = sqrt(eps(class(mu)));
    alpha = (1-f)*calphamax(index)+f*calphamin(index); % farthest edge, avoiding roundoff error
  end
  if isnan(alpha)
    warning('cn:linesearch','alpha is NaN')
  end
end
