function [pars] = sine_fit(means,raw)
% SINE_FIT fits a sinusoid function to data using 'fmincon', with parameter bounds
%	It calls sine_err.m to comput the error for each set of params.  THe function 
%   evaluated is given by sine_func.m

global Data RawData;

%allows calculation of error for both raw data and mean data
Data = means;
RawData = raw;

% first, generate some initial parameter guesses
[max_val max_indx] = max(Data(:,2));
[min_val min_indx] = min(Data(:,2));
[max_x max_x_indx] = max(Data(:,1));
[min_x min_x_indx] = min(Data(:,1));
N_values = length(Data(:,1));

q(1) = mean(Data(:,2));
q(2) = 0.5*(max_val - min_val);
q(3) = 1.0/(max_x - min_x); %one cycle per x-range
q(4) = 0.5;

%get a good estimate of the starting frequency from the Fourier Transform
FFT_LENGTH = 128;
FT = abs(fft(Data(:,2) - mean(Data(:,2)), FFT_LENGTH));
FT2 = FT(1:length(FT)/2);
[max_ampl max_indx] = max(FT2);
%compute the disparity interval
dx = (max_x - min_x)/(N_values-1);
freq = (max_indx/length(FT2))/(2*dx);
q(3) = freq;

% do a quick search for a good starting phase
errors = [];
N=20;
for i=1:N
   q(4) = -pi + 2*pi*(i-1)/N;  
   errors(i) = sine_err(q);
end
[min_err min_indx] = min(errors);
q(4) = -pi + 2*pi*(min_indx-1)/N;


A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=[0; 0; 0.1/(max_x - min_x); -2*pi];  %lower bounds
UB=[1.5*max_val; 1.35*(max_val - min_val); 5.0/(max_x - min_x); 2*pi]; %upper bounds

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');
%OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000);

N_reps = 40;
wiggle = 0.2;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    testpars{j} = fmincon('sine_err',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = sine_err(testpars{j});
end
%err

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};
return;