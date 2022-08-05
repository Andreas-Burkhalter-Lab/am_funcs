function y = deconv_exponential(y,timeconst,n0)
% DECONV_EXPONENTIAL: estimate underlying signal from filtered, noisy version
%
% Suppose we have a signal y. Filter this signal with an exponential
% filter, and then add some noise. This filtered, noisy signal will be
% called yn. The goal is to attempt to reconstruct y given yn, the
% exponential decay constant, and the noise amplitude. Our estimate of y is
% called yd, the "deconvolved" yn.
%
% Syntax:
%   yd = deconv_exponential(yn,timeconst,n0)
% where
%   yn is the signal (a vector or a matrix; each column is one signal);
%   timeconst is the time constant of the exponential decay, in samples;
%   n0 is the noise level, in units where 1 is the magnitude of the signal;
% and
%   yd is the estimated underlying signal that produced y.
%
% Example:
% % First, make a filtered, noisy signal
% y0 = zeros(1,1000); y0(100) = 1; y0(200) = 1; %(the underlying signal)
% alpha = 1/100; a(2) = -exp(-alpha); a(1) = 1; b = 1; y = filter(b,a,y0);
% yn = y + 0.1*randn(size(y)) - (1:1000)/1000;  %(add some photobleaching)
%
% % Now do the deconvolution
% yd = deconv_exponential(yn,1/alpha,0.03);
%
% % Plot
% figure; plot(yn)
% figure; plot(yd)



% Copyright 2006 by Timothy E. Holy

  [n,k] = size(y);
  if (n == 1)
    y = y(:);
    [n,k] = size(y);
  end
  % Calculate filtering coefficients
  alpha = 1/timeconst;
  beta = alpha*sqrt(1+1/n0^2);
  a = [1 -exp(-beta)];
  b = 1;
  
  % Calculate filter initial conditions
  nfilt = max(length(a),length(b));
  a(end+1:nfilt) = 0;
  b(end+1:nfilt) = 0;
  rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
  cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
  data = [1+a(2) a(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
  sp = sparse(rows,cols,data);
  zi = sp \ ( b(2:nfilt).' - a(2:nfilt).'*b(1) );

  % Pad the beginning and end
  nfact = 3*(nfilt-1);
  ypad = [2*y(ones(1,nfact),:)-y(nfact+1:-1:2,:); y; 2*y(n*ones(1,nfact),:)-y(n-1:-1:n-nfact,:)];

  % Filter
  yff = filter(b,a,ypad,zi*ypad(1,:));
  yfb = filter(b,a,ypad(end:-1:1,:),zi*ypad(end,:));
  yfb = yfb(end:-1:1,:);
  
  % Remove the padding
  yff = yff(nfact+1:end-nfact,:);
  yfb = yfb(nfact+1:end-nfact,:);
  
  % Compute the linear coefficients
  coeff = 1/(1+n0^2) - 1/(n0*sqrt(1+n0^2));
  coefb = 1/(1+n0^2) + 1/(n0*sqrt(1+n0^2));
  
  % Finally, compute the filtered product
  y = coeff*yff + coefb*yfb;
  