function [crange,trs] = resample_gradient(gc,N)
% resample_gradient: set resampling times to yield a common linear gradient
%
% Given a chromatography gradient schedule of the following form:
%    0  c0
%    t1 c1
%    t2 c2
%    ...
% where t_i is the time at which the concentration of B becomes c_i, the
% task is to "standardize" all samples, even if they employed a different
% schedule.
%
% This function finds the monotonic-increasing portion of each of a list of
% gradients, and computes a set of times so that each "run", if resampled
% at those times, would correspond to a linearly-increasing concentration
% axis.
%
% Syntax:
%   [crange,trs] = resample_gradient(gc,N)
% where
%   gc is a cell array with n_files elements, each a gradient table
%     arranged as shown above;
%   N is the # of samples one wishes to have in the output
% and
%   crange is [cmin cmax], the range of concentrations common to all of the
%     files;
%   trs is n_files-by-N, each row corresponding to a file. It contains the
%     list of times at which the concentration was equal to
%     linspace(cmin,cmax,N).

% Copyright 2010 by Timothy E. Holy

  n_files = length(gc);
  % For each file, find the beginning and end of the monotonic rise
  gsnipc = cell(1,n_files);
  indexsnip = zeros(2,n_files);
  for k = 1:n_files
    g = gc{k};
    dg = diff(g,1,1);
    indxstart = find(dg(:,2) > 0,1,'first');
    if isempty(indxstart)
      error('There is no gradient for this file');
    end
    n = find(dg(indxstart:end,2) < 0,1,'first');
    if isempty(n)
      indxend = size(dg,1);
    else
      indxend = n+indxstart-1;
    end
    while(indxend > indxstart && g(indxend,2) == g(indxend-1,2))
      indxend = indxend-1;
    end
    g = g(indxstart:indxend,:);
    gsnipc{k} = g;
    indexsnip(:,k) = [indxstart; indxend];
  end
  % Within the monotonic period, find range of concentration that is common
  % to all files
  cmin = max(cellfun(@(e) min(e(:,2)),gsnipc));
  cmax = min(cellfun(@(e) max(e(:,2)),gsnipc));
  crange = [cmin cmax];
  % Find sampling times that correspond to a linear mapping of this range
  trs = zeros(n_files,N);
  ctarget = linspace(cmin,cmax,N);
  for k = 1:n_files
    trs(k,:) = interp1(gsnipc{k}(:,2),gsnipc{k}(:,1),ctarget);
  end
