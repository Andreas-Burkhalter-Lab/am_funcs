function [pars,selected_freq, Accounted_Var] = gaborgauss_sumfit(means,raw,fixed_param_flags,fixed_param_values, CCG_std)
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
global CCG_error
%allows calculation of error for both raw data and mean data
Data = means;
RawData = raw;
CCG_error = CCG_std;
% first, generate some initial parameter guesses
[max_val	max_indx] = max(Data(:,2));
[min_val	min_indx] = min(Data(:,2));
[max_time max_time_indx] = max(Data(:,1));
[min_time min_time_indx] = min(Data(:,1));
N_values = length(Data(:,1));

q(1) = mean(Data(:,2));         %offset
q(2) = 0.5*(max_val - min_val); %amplitude of gabor
q(3) = 0.0;                     %center
q(4) = 1;                     %size of gabor
q(5) = 0.5;                     %freq
q(6) = 0.5*(max_val - min_val); %amplitude of gaussian
q(7) = 2;                       %exponent of gabor
q(8) = 2;                      % size of gaussian

%get a good estimate of the starting frequency from the Fourier Transform
[freq, ampl] = FourierTransform_1D( 1/1000 *[1 : (length(Data(:,2)) - 1)], Data(2:end, 2), length(Data(:,2) ) - 1, 1, 0);  
[max_ampl max_indx] = max(ampl(freq > 20) );
selected_freq = freq(freq > 20);
selected_freq = selected_freq(max_indx);
if (selected_freq > 100)
    selected_freq = 50;
end

q(5) = selected_freq;     
A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=[0;0;-10;10;10; 0; 2; 0.1];
UB=[1.5*max_val;2*(max_val - min_val);10;2*(max_time - min_time);100; 1.5*(max_val - min_val); 5.0; 20];

%Here we set the starting parameter values and bounds to the fixed values, if requested by the flags
q(logical(fixed_param_flags)) = fixed_param_values(logical(fixed_param_flags));
LB(logical(fixed_param_flags)) = fixed_param_values(logical(fixed_param_flags));
UB(logical(fixed_param_flags)) = fixed_param_values(logical(fixed_param_flags));
%[q' LB UB]

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');

N_reps = 20;
wiggle = 0.3;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    rand_factor(logical(fixed_param_flags)) = 1.0;  %don't randomize for parameters that are fixed
    temp_q = q' .* rand_factor;
    temp_q(5) = q(5); %don't override the frequency estimate
    testpars{j} = fmincon('gaborgauss_sumerr',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = gaborgauss_sumerr(testpars{j});
end
%err

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};

% fit_cross_correlogram = real (gaborgauss_sumfunc(Data(:,1), pars) ) ;
CCG_Var = sum ((  (Data(:,2) - mean(Data(:,2)) )./ CCG_error'  ).^2    );
% GaborFit_Var= sum((fit_cross_correlogram - Data(:,2)   ).^2 );
% Accounted_Var = 100*( 1 - (GaborFit_Var /CCG_Var ) );
% if Accounted_Var < 15
%     %If fit curve accounts for only 15% of variance, set fit to flat
%     %line
%     %keep baseline the same, just use flat function
%     pars(2:8) = 0;
%     pars(1) = mean((Data(:,2)));
% end

Accounted_Var = 100*( 1 - (min_err/CCG_Var ) );
if  (Accounted_Var < 15)
%     %If fit curve accounts for only 15% of variance, set fit to flat
%     %line
%     %keep baseline the same, just use flat function
    pars(2:8) = 0;
    pars(1) = mean((Data(:,2)));
end



return;