function pars = sinfit(means,raw)

global Data RawData;

Data = means;
RawData = raw;

% first, generate some initial parameter guesses
[max_val	max_indx] = max(Data(:,2));
[min_val	min_indx] = min(Data(:,2));
[max_disp max_disp_indx] = max(Data(:,1));
[min_disp min_disp_indx] = min(Data(:,1));N_values = length(Data(:,1));


%initialize parameters--------------------------------------------------
%amplitude
q(1) = (max_val-min_val)/2;

%frequency
q(2) = 1;

%baseline
q(4) = mean(Data(:,2));

% do a quick search for a good starting phase
errors = [];
N=20;
for i=1:N
   q(3) = -pi + 2*pi*(i-1)/N;  
   errors(i) = sinerr(q);
end
[min_err min_indx] = min(errors);
q(3) = -pi + 2*pi*(min_indx-1)/N;

%initialize parameters--------------------------------------------------

A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=[ 0; 1; -2*pi; 0];
UB=[ 500; 1; 2*pi; 500];

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'final', 'MaxFunEvals', 1000);

N_reps = 40;
wiggle = 0.2;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    temp_q(2) = q(2); %don't override the frequency estimate
    testpars{j} = fmincon('sinerr',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = sinerr(testpars{j});
 end

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};
return;