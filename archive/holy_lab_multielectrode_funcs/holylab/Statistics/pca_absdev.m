function [C,R,iter] = pca_absdev(X,C0)
% PCA_ABSDEV: a variant of PCA robust to outliers
% Syntax:
%   C = pca_absdev(X,n_components)
%   C = pca_absdev(X,C0)
%   [C,R] = pca_absdev(...)
%   [C,R,iter] = pca_absdev(...)
% where
%   X is a npoints-by-d matrix
%

  [N,d] = size(X);
  if (isscalar(C0) && d > 1)
    n_components = C0;
    % Build up a set of initial directions by choosing random elements
    % from X, until we've got a subspace of dimension n_components
    % (this is espec. necessary if X sometimes comes from bootstrap
    % resampling, so we can't just assume that n_components different
    % random draws from X will give us what we need)
    [C0,R] = pca_absdev_initialize(X,[],[],n_components);
    while (size(C0,2) < n_components)
      [C0,R] = pca_absdev_initialize(X,C0,R,n_components);
    end
  end

    randIndex = round(N*rand(1,n_components)+0.5);
    C0 = X(randIndex,:)';
    % Now make sure we ended up with an orthonormal subset
    [C0,R] = qr(C0);
    diagR = abs(diag(R));
    goodIndex = find(diagR > sqrt(eps)*diagR(1));
    C0 = C0(:,goodIndex);
    while (size(C0,2) < n_components)
      
  end
  [C0,tmp] = qr(C0,0);
  [C,width] = pca_absdev_iter(X,C0);
  iter = 1;
  Cdiff = max(max(abs(C-C0)));
  while (Cdiff > eps)
    fprintf('max change: %g',max(max(abs(C-C0))))
    C0 = C;
    [C,width] = pca_absdev_iter(X,C0);
    iter = iter+1;
    Cdiffold = Cdiff;
    Cdiff = max(max(abs(C-C0)));
    fprintf(' d(max change): %g\n',Cdiff-Cdiffold)
    %if (abs(Cdiff - Cdiffold) < eps)
    %  fprintf('Triggered average!\n')
    %  C = (C+C0)/2;
    %end
  end
  width = width/N;
  width'

function [C,width] = pca_absdev_iter(X,C0)
% supply X as a N-by-d matrix (observations are rows), C0 as a column
% matrix
  [N,dim] = size(X);
  proj = X*C0;
  S = sign(proj);
%  w = sum(abs(proj));
  P = X'*S;
  L0 = sum(sum(abs(proj)));
%  Pw = repmat(w,dim,1) .* P;
  %[C,lambda] = poldec(P);
  [C,R] = qr(P,0);
  %lambda = sqrtm(P'*P);
  %C = P/lambda;
  width = diag(R);
  %L = trace(Z*V');
  L = sum(width);
  %if (L < L0)
  %  fprintf(' reflect ');
  %  C = -C0;
  %end
  fprintf(' L0 %g, L %g ',L0,L);
  %if (L < L0)
  %  keyboard
  %end
