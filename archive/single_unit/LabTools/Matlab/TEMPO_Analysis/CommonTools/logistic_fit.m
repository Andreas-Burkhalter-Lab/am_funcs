function [alpha, beta] = logistic_fit(data)
% LOGISTIC_FIT fits a logistic function to data using maximum likelihood
%	maximization under binomial assumptions(?)  It uses LOGISTIC_FUNC for
%	error calculation.  Data must be in 3 columns: x, %-correct, Nobs
%	usage: [alpha, beta] = logistic_fit(data)
%		alpha and beta are the threshold and slope parameters.
global Data;
Data = data;

% generate guess
q = ones(2,1);

% use linear interpolation to guess alpha
% next line takes x and %-cor columns, flips them and interpolates along
% x to find the value corresponding to .8 correct.  The interpolation
% requires monotonic %-cor column, so we sort the matrix 1st.
% Remember, it's just a guess.
%a = find(data(:,2)>.65 & data(:,2)<.95);
%if isempty(a)
%  q(1,1) = mean(data(:,1));
%else
%  q(1,1) = mean( data(a,1) );
%end  

%get a starting threshold estimate by testing a bunch of values
%this should work better than above.
N = 30;
for i=1:N
   q(1,1) = (max(data(:,1))/N)*i;
   errors(i) = logistic_func(q);
end
[min_err min_indx] = min(errors);
q(1,1) = (max(data(:,1))/N)*min_indx;

OPTIONS = OPTIMSET('MaxIter', 10000,'MaxFunEvals',2000);
quick = fminsearch('logistic_func',q);
% quick
alpha = quick(1,1);
beta = quick(2,1);


