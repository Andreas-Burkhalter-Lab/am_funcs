function [pars] = diff_erf_fit(means,raw)
% diff_erf_fit fits a Difference-of-Error function to data using 'fmincon', with parameter bounds
%	It calls diff_erf_err.m to compute the error for each set of params.  THe function 
%   evaluated is given by diff_erf_func.m

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

q(1) = max_val; % amplitude of excitatory erf
q(2) = 3;   %size constant of excitatory erf
q(3) = max_val - Data(N_values,2); %amplitude of inhibitory erf
q(4) = 3; %size constant of inhibitory erf (above and beyond excitatory erf)
q(5) = min_val; %baseline firing rate

A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=[0; 0.01; 0; 0.001; -max_val/10];  %lower bounds
UB=[5*max_val; 100; 5*max_val; 100; max_val]; %upper bounds

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');
%OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000);

N_reps = 10;
wiggle = 0.2;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    testpars{j} = fmincon('diff_erf_err',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = diff_erf_err(testpars{j});
end
%err

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};
return;