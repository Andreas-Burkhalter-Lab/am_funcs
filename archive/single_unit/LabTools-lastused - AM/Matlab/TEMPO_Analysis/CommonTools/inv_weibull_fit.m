function [alpha, beta] = inv_weibull_fit(data)
% WEIBULL_FIT fits a weibull function to data using maximum likelihood
%	maximization under binomial assumptions.  It uses WEIBULL_FUNC for
%	error calculation.  Data must be in 3 columns: x, %-correct, Nobs
%	usage: [alpha, beta] = weibull_fit(data)
%		alpha and beta are the threshold and slope parameters.
global Data;
Data = data;
% generate guess
% slopes of around 4 common for fit curves
q = 4 * ones(2,1);

% use linear interpolation to guess alpha
% next line takes x and %-cor columns, flips them and interpolates along
% x to find the value corresponding to .8 correct.  The interpolation
% requires monotonic %-cor column, so we sort the matrix 1st.
% Remember, it's just a guess.
a = find(data(:,2)>.7 & data(:,2)<.9);
if isempty(a)
  q(1,1) = mean(data(:,1));
else
  q(1,1) = mean(a);
end  
tmpdata = data;
% kludge:
% table1 barfs if the column of numbers it is handed
% contains duplicates (complains about the numbers
% being "non-monotonic") so we add a small random
% number to duplicates.  Blech!!
for i = 1:size(data,1)
  mask = logical(ones(size(data,1),1));
  mask(i) = 0;
  if any(data(i,2) == data(mask,2))
    tmpdata(i,2) = tmpdata(i,2) + rand/10000;
  end
end  



if (.7 > max(data(:,2)))
  q(1,1) = max(data(:,1));
else
  q(1,1) = table1(fliplr(sort(tmpdata(:,1:2))),.7);
end
trace = 0;
tol = .0001;
quick = fmins('inv_weibull_func',q,[trace tol]);
% quick
alpha = quick(1,1);
beta = quick(2,1);


