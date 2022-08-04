function [pars,freq] = gaborfit_for_SurroundModel(means,raw,fixed_param_flags,fixed_param_values, uncorr)
% GABORFIT fits a Gabor function to data using 'fmincon', with parameter bounds
%	It calls gaborerr.m to comput the error for each set of params.  THe function 
%   evaluated is given by gaborfunc.m
% 'means' is an Nx2 array, where the first column is each unique disparity and the second column is the mean response at that disparity
% 'raw' is an Mx2 array where the first column is the disparity for each trial and the second column is the response for each trial
% Hence, either the mean responses or the raw trial responses can be fit depending on what is selected (by hard coding) in gaborerr.m
% 'fixed_param_flags' should contain zeros for each function parameter that you want to be varied (fill with zeros to vary all params)
% 'fixed_param_values' should contain the desired values for any of the parameters that you want to be held fixed during fitting
%
% this function returns the best-fitting parameters, as well as the disparity frequency derived from the FT of the raw data

global Data RawData;

%allows calculation of error for both raw data and mean data
Data = means;
RawData = raw;

% first, generate some initial parameter guesses
[max_val	max_indx] = max(Data(:,2));
[min_val	min_indx] = min(Data(:,2));
[max_disp max_disp_indx] = max(Data(:,1));
[min_disp min_disp_indx] = min(Data(:,1));
N_values = length(Data(:,1));

q(1) = (max_val-min_val)/2;
q(2) = 0.5*(max_val - min_val);
q(3) = 0.0;
q(4) = 0.5; % 
q(5) = 0.5; % 
q(6) = -1.5;

BR_LB = min_val - (.2*min_val);
BR_UB = max_val + (.2*max_val);


%get a good estimate of the starting frequency from the Fourier Transform
FFT_LENGTH = 512;
FT = abs(fft(Data(:,2) - mean(Data(:,2)), FFT_LENGTH));
FT2 = FT(1:length(FT)/2);
[max_ampl max_indx] = max(FT2);
%compute the disparity interval
dx = (max_disp - min_disp)/(N_values-1);
freq = (max_indx/length(FT2))/(2*dx);
q(5) = freq;

% do a quick search for a good starting phase
errors = [];
N=20;
for i=1:N
   q(6) = -pi + 2*pi*(i-1)/N;  
   errors(i) = gaborerr(q);
end
[min_err min_indx] = min(errors);
q(6) = -pi + 2*pi*(min_indx-1)/N;

A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=[BR_LB;0;min_disp;0.1;0.9*freq;-2*pi];
UB=[BR_UB;1.35*(max_val - min_val);max_disp;(max_disp - min_disp)/2;1.1*freq;2*pi];

%Here we set the starting parameter values and bounds to the fixed values, if requested by the flags
q(logical(fixed_param_flags)) = fixed_param_values(logical(fixed_param_flags));
LB(logical(fixed_param_flags)) = fixed_param_values(logical(fixed_param_flags));
UB(logical(fixed_param_flags)) = fixed_param_values(logical(fixed_param_flags));
%[q' LB UB]

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');
%OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000);

N_reps = 40;
wiggle = 0.3;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    rand_factor(logical(fixed_param_flags)) = 1.0;  %don't randomize for parameters that are fixed
    temp_q = q' .* rand_factor;
    temp_q(5) = q(5); %don't override the frequency estimate
    testpars{j} = fmincon('gaborerr',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = gaborerr(testpars{j});
end
%err

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};
return;