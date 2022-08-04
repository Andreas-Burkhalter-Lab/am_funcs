function [pars] = hyperbolicratio_fit(means,raw);
% Fits contrast data with a hyperbolic ratio function.  Based on gaussfit
% protocol.  CMA 05/24/06
%         (Rmax * c^n)/
% R(c) = (sigma^n + c^n)
% Where c is stimulus contrast, R is the cell's response, and Rmax, sigma,
% and n are free parameters (SEE:  "A Tonic Hyperpolarization Underlying 
% Contrast Adaptation in Cat Visual Cortex" SCIENCE Vol 276, 9 May 1997)

%%% 3/14/16 AM edited to make compatible with newer matlab versions 

global Data RawData;

N=30;

%allows calculation of error for both raw data and mean data
Data = means;
RawData = raw;

% first, generate some initial parameter guesses
[max_val max_indx] = max(Data(:,2));
[min_val min_indx] = min(Data(:,2));
[max_x max_x_indx] = max(Data(:,1));
[min_x min_x_indx] = min(Data(:,1));
N_values = length(Data(:,1));

q(1) = max_val;             % Rmax
q(2) = 40;                  % sigma
q(3) = 2;                   % n
q(4) = min_val;             %DC offset

A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=[0; 0; 0; 0];  %lower bounds
UB=[1.0*max_val; 100; 10; max_val]; %upper bounds

%%OPTIONS = OPTIMSET('fmincon');
%%OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');
%OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000);
OPTIONS = optimset('MaxIter', 5000); % 3/14/16 AM replaced preceding 2 lines with this line

N_reps = 40;
wiggle = 0.2;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    testpars{j} = fmincon('hyperbolicratio_error', temp_q, A, B, Aeq, Beq, LB, UB, NONLCON, OPTIONS);
    err(j) = hyperbolicratio_error(testpars{j});
end
%err

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};

return;