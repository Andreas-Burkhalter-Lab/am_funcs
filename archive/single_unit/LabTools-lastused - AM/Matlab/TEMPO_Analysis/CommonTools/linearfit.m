function [pars] = linearfit(means,fixed_param_flags,fixed_param_values)
% LINEARFIT fits a linear function to data using 'fmincon', with parameter bounds
%	It calls linearerr.m to compute the error for each set of params.  The function 
%   evaluated is given by linearfunc.m 
%   Use this in cases where you want to constrain parameters, otherwise you can use REGRESS. TU 01/03/02

global Data;

%allows calculation of error for both raw data and mean data
Data = means;

% first, generate some initial parameter guesses
max_y = max(Data(:,2));
min_y = min(Data(:,2));
max_x = max(Data(:,1));
min_x = min(Data(:,1));

q(2) = (max_y - min_y) / (max_x - min_x);
q(1) = Data(1,2) - q(2) * Data(1,2);

A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[]; 
LB=[-999999 , -999999 ]; 
UB=[999999, 999999];

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
    testpars{j} = fmincon('linearerr',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = linearerr(testpars{j});
end
%err

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};
return;