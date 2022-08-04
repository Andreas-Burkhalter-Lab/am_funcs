function [pars] = gammafit(means,raw)
% GAMMAFIT fits a Gaussian function to data using Simplex algorithm.
%	It calls gaussfunc.m to comput the error for each set of params
%	Params are as follows: Ro (q1) =base rate, K (q2) =amplitude, x0 (q3) =center, a (q4) =size
global Data RawData;

%allows calculation of error for both raw data and mean data
Data = means;
RawData = raw;

% first, generate some initial parameter guesses
[max_val	max_indx] = max(Data(:,2));
[min_val	min_indx] = min(Data(:,2));
[max_x max_x_indx] = max(Data(:,1));
[min_x min_x_indx] = min(Data(:,1));
N_values = length(Data(:,1));

%these starting conditions generally work well
q(1) = min_val;
q(2) = max_val - min_val;
q(3) = .1; 
q(4) = 0;  
q(5) = .5; 

%do a coarse search to find good starting values for parameters 3-5
min_err = gamma_err(q);
for i = 0:0.05:0.8
    for j = -4:0.1:2
        for k = 0:0.2:5
            qtemp = q;
            qtemp(3) = i;
            qtemp(4) = j;
            qtemp(5) = k;
            temp_err = gamma_err(qtemp);
            if (temp_err < min_err)
                min_err = temp_err;
                q(3) = i;
                q(4) = j;
                q(5) = k;
            end
        end
    end
end
            
A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=[0;0;0;-20;0];
UB=[1.5*min_val;1.5*(max_val - min_val);50;20;100];

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');
%OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000);

N_reps = 10;
wiggle = 0.3;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    testpars{j} = fmincon('gamma_err',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = gamma_err(testpars{j});
end
%err

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};
return;
