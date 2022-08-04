function [pars, fit_err] = gammafit2(data, q)
% GAMMAFIT fits a Gaussian function to data using Simplex algorithm.
%	It calls gaussfunc.m to comput the error for each set of params
%	Params are as follows: Ro (q1) =base rate, K (q2) =amplitude, x0 (q3) =center, a (q4) =size
global Data;
Data = data;

% first, generate some initial parameter guesses, which is pretty easy for a Gaussian
[max_val	max_indx] = max(Data(:,2));
[min_val	min_indx] = min(Data(:,2));

%These special starting conditions work well for curves that peak out at slow speeds
%and have higher responses at faster speeds
%q(1) = 226;
%q(2) = -107;
%q(3) = .28;
%q(4) = .014;
%q(5) = 2.12;

%these starting conditions generally work well
%q(1) = min_val;
%q(2) = max_val - min_val;
%q(3) = 5; %log(max_val - Data(1,2))/log(Data(max_indx,1) - Data(1,1)); % exponential decay
%q(4) = 0;  % check this offset need a better guess
%q(5) = 5; %log(max_val - Data(1,2))/log(Data(max_indx,1) - Data(1,1)) % exponential rise
start(3) = q(3);
start(4) = q(4);
start(5) = q(5);
rand('state',sum(100*clock));

OPTIONS(14) = 1400;
fit_err = 0;
iter = 1;
for iter = 1:10  % must revise still have problems
   quick(iter,:) = fmins('gammafunc',q, OPTIONS);
   err(iter) = gammafunc(quick(iter,:));	  
   q(3) = start(3) + 10*randn;   
   q(5) = start(5) + 10*randn; 
end
err(err == 0) = NaN; 			% do not consider meaningless fits, set error of 0 to NaN
[fit_err fit_i] = min(err);
pars = quick(fit_i,:);
