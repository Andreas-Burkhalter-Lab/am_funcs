function [pars] = gauss2Dfit(means,raw)
% gaussfit fits a sinusoid function to data using 'fmincon', with parameter bounds
%	It calls gausserr.m to comput the error for each set of params.  THe function 
%   evaluated is given by gaussfunc.m

global Data RawData;

%allows calculation of error for both raw data and mean data
Data = means; %[px py spikes]
RawData = raw; %[plot_x' plot_y' spikes]

% first, generate some initial parameter guesses
[max_val max_indx] = max(Data(:,3));
[min_val min_indx] = min(Data(:,3));
[max_y max_y_indx] = max(Data(:,2));
[min_y min_y_indx] = min(Data(:,2));
[max_x max_x_indx] = max(Data(:,1));
[min_x min_x_indx] = min(Data(:,1));
N_values = length(Data(:,1));

q(1) = min_val;
q(2) = (max_val - min_val);
q(3) = Data(max_indx,1);
q(4) = 0.2*(max_x - min_x);
q(5) = Data(max_indx,1);
q(6) = 0.2*(max_y - min_y);


A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=[0; -1.5*(max_val - min_val); min_x; 0.05*(max_x - min_x); min_y; 0.05*(max_y - min_y)];  %lower bounds
UB=[1.35*max_val; 1.5*(max_val - min_val); max_x; 2*(max_x - min_x); max_y; 2*(max_y - min_y)]; %upper bounds

OPTIONS = OPTIMSET('fmincon'); %
%OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000);

N_reps = 20;
wiggle = 0.2;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    testpars{j} = fmincon('gauss2Derr',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = gauss2Derr(testpars{j});
end
%err

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};
return;