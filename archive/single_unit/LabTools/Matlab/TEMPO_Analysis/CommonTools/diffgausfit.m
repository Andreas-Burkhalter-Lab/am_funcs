function [pars] = diffgausfit(data)
%DIFFGAUSFIT - performs a difference of gaussian fit.  this is to aid in creating a fit of data to receptive
%fields with surrounds

global Data
Data = data;

% first, generate some initial parameter guesses, which is pretty easy for a Gaussian
[max_val	max_indx] = max(Data(:,2));
[min_val	min_indx] = min(Data(:,2));

[max_size max_indx] = max(Data(:,1));

start(1) = max_val;
start(2) = 2;
start(3) = .75*start(1);
start(4) = 6;

q(1) = max_val;
q(2) = 2;
q(3) = .75*q(1);
q(4) = 6;

options(14) = 200*4;

rand('state',sum(100*clock));

for iter = 1:10  % code taken from ben's gammafit function
   quick(iter,:) = fmins('diffgausfunc',q, options);
   err(iter) = diffgausfunc(quick(iter,:));	  
   q(1) = start(1) + 10*randn;   
   q(2) = start(2) + 10*randn; 
   q(3) = start(3) + 10*randn;   
   q(4) = start(4) + 10*randn; 
end
err(err == 0) = NaN; 			% do not consider meaningless fits, set error of 0 to NaN
err
quick
[fit_err fit_i] = min(err)
pars = quick(fit_i,:);

