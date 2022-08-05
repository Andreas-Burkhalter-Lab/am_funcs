function [imz,it] = msprof_factor(ic)
% MSPROF_FACTOR: separate MS sweeps into products of time & m/z
% The ion current ic(mz,t) is a function of m/z and time t, where time
% enters because the sample is injected and moves through the tubing
% before hitting the detector.  We'd like to separate out the
% "background" or "contaminant" component from the "injected" component,
% and eliminate the time-dependence of the signal.
%
% The idea here is to write the ion current in the following way:
%     ic(mz,t) = sum_j imz_j(mz) * it_j(t),
% i.e., to seek a representation as a sum of separable functions.  Both
% imz_j and it_j should be non-negative, as negative material has no
% meaning. So this is like nonnegative matrix factorization, except
% we will restrict the "background" to have constant time-dependence.
%
% Syntax:
%   [imz,it] = msfactor(ic)
% where
%   ic is a n_mz-by-n_sweeps matrix of the ion currents recorded as a
%     function of m/z and time;
% and
%   imz,it are non-negative matrices so that 
%      ic "=" imz*it.
% imz is n_mz-by-2, it is 2-by-n_sweeps.  ("=" means approximately equal)
% The first column of imz (row of it) refers to the background, and row 1
% of it is constant. The 2nd column of imz is the spectrum of the injected
% material, and the 2nd row if it is its timecourse.
%
% The timecourses are normalized to have sum = 1, so that the spectrum is
% proportional to the total quantity of injected material.
%
% While we pursue our own (better, for this 2-component problem) algorithm
% here, see paper on non-negative matrix factorization (NMF) by Seung & Lee.
%
% See also: MSPROF_PARSE.
  
% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'force_const')
    options.force_const = 0;
  end
  % Temporarily eliminate rows with all zeros
  iz = sum(ic,2) == 0;
  icnz = ic(~iz,:);
  % Factor using non-negative matrix factorization, forcing the time for
  % the first component to be constant
  [m,n] = size(icnz);
  r = 2;
  imz=rand(m,r);
  it=rand(r,n);
  % Set up the stuff needed for iteration
  errold = sum(sum((icnz - imz*it).^2));
  isdone = 0;
  iter = 0;
  fprintf('%d: %g\n',iter,errold);
  while ~isdone
    % Solve for the times. We can do this directly because there's only
    % one unknown row to the time matrix; the top row is 1/n.
    num = sum((icnz - repmat(imz(:,1)/n,1,n)).*repmat(imz(:,2),1,n));
    denom = sum(imz(:,2).^2);
    t = max(0,num/denom); % Non-negative. This isn't a hack, it solves
                          % LSQ exactly in this case.
    t = t/sum(t);         % Normalize
    it = [1/n + zeros(1,n); t];
    % Solve for the m/z spectra.  We'll also do this by hand. It's a set
    % of uncoupled equations in 2 variables, so we can calculate the
    % solution explicitly.
    % First calculate the solution in the absence of constraints
    b = it*icnz';
    H = it*it';  % Hessian
    imzopt = H\b;
    % Find the ones that have at least one of the two coordinates in the
    % forbidden region, and then restrict that coordinate to being on the
    % axis line
    if any((imzopt(1,:) < 0) & (imzopt(2,:) < 0))
      % This means I don't understand what's going on
      error('Oops, confusing')
    end
    indxbad1 = find(imzopt(1,:) < 0);
    imzopt(1,indxbad1) = 0;
    imzopt(2,indxbad1) = b(2,indxbad1)/H(2,2);
    indxbad2 = find(imzopt(2,:) < 0);
    imzopt(2,indxbad2) = 0;
    imzopt(1,indxbad2) = b(1,indxbad2)/H(1,1);
    imz = imzopt';
    err = sum(sum((icnz - imz*it).^2));
    if (err > (1-1e-5)*errold)
      isdone = 1;
    else
      errold = err;
      iter = iter+1;
      fprintf('%d: %g\n',iter,err);
    end
  end
  % Restore original size
  tmp_imz = imz;
  imz = zeros(size(ic,1),r);
  imz(~iz,:) = tmp_imz;