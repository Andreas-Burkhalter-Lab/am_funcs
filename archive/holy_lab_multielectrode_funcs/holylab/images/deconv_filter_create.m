function f = deconv_filter_create(m,noise)
% DECONV_FILTER_CREATE: create a deconvolution filter
% This calculates a Wiener filter for one-dimensional data, in the
% real-time domain. Thus, it is not well-suited for very large problems.
%
% Syntax:
%   f = deconv_filter_create(m,noise)
% where
%   m contains the measured pixel intensities (must be odd in length);
%   noise is the combination of shot noise + CCD readout noise, in
%     digital numbers (default 0);
% and
%   f is the appropriate deconvolution filter with respect to digital
%     numbers.

% Copyright 2006 by Timothy E. Holy

  % Check to make sure that the length of m is odd
  n = length(m);
  halfn = round((n-1)/2);
  if (2*halfn+1 ~= n)
    error('The length of m must be odd');
  end
  
  if (nargin < 2)
    noise = 0;
  end
  
  mc = conv(m,m);
  summ = sum(m);
  m_rhs = summ*m(end:-1:1)';
  for i = 1:length(m_rhs)
    for j = 1:length(m_rhs)
      A(i,j) = mc(n + i - j);
    end
  end
  % Add the noise terms
  for i = 1:length(m_rhs)
    A(i,i) = A(i,i) + noise^2;
  end

  f = (A\m_rhs)';

  