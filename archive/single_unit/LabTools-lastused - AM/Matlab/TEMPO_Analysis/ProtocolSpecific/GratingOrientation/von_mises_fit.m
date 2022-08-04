function [pars] = von_mises_fit(means,raw)
%Fits orientation data with a von_mises function - CMA 10/01/06
%M(x) = A1*exp{k1[cos2*(theta-psi1)-1]} + A2*exp{k2[cos2*(theta-psi2)-1]} + DC;
%(sum of two von mises fuctions)
%A is the value of the function at the preferred orientation, psi, and k is
%a width parameter
%
%q(1) is A1, q(2) is k1, q(3) is psi, q(4) is A2, q(5) is k2, q(6) is psi2,
%q(7) is a DC offset term
%SEE: (Swindale NV: Biological Cybernetics 78:45-56, 1998)
%%% AM edited 10/2/15 to fix for new version of Matlab/Windows

global Data RawData;

%allows calculation of error for both raw data and mean data
Data = means;
RawData = raw;

%let the center of the first peak be at the location of the maximum value
[m n] = size(Data);
[max_val max_indx] = max(Data(:,2));
[min_val min_indx] = min(Data(:,2));

%try to find the second peak of the data under the constraint that the
%second peak can't be within 30 degrees of the first peak
if max_indx == 1
    NewData = Data;
    NewData(1:2,2) = 0;
    NewData(m,2) = 0;
elseif max_indx == n
    NewData = Data;
    NewData((m-1):n,2) = 0;
    NewData(1,2) = 0;
else
    NewData = Data;
    NewData((max_indx-1):(max_indx+1),2) = 0;
end

[max_val2 max_indx2] = max(NewData(:,2));
    

q(1) = max_val-min_val;             % A1
q(2) = 10;                  % k1
q(3) = Data(max_indx,1);    % psi1
q(4) = max_val2-min_val;    % A2
q(5) = 10;                  % k2
q(6) = q(3)+180;
q(7) = min_val;             % DC


A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=[0.4*(max_val-min_val); 0; q(3)-15; 0.4*(max_val-min_val); 0; q(6)-15; 0];  %lower bounds
UB=[1.2*max_val; 100; q(3)+15; 1.2*max_val; 100; q(6)+15; min_val]; %upper bounds

% OPTIONS = optimset('fmincon'); %%%% AM modified 'OPTIMSET' to 'optimset' 10/1/15
% OPTIONS = optimset('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off'); %%%% AM modified 'OPTIMSET' to 'optimset' 10/1/15
OPTIONS = optimset('MaxIter',5000); %% AM 10/2/15 replaced above two lines with this line as per Aki's instructions
%OPTIONS = optimset('LargeScale', 'off', 'LevenbergMarquardt', 'on',   %%%% AM modified 'OPTIMSET' to 'optimset' 10/1/15
%'MaxIter', 5000); 

N_reps = 40;
wiggle = 0.2;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    testpars{j} = fmincon('von_mises_error', temp_q, A, B, Aeq, Beq, LB, UB, NONLCON, OPTIONS);
    err(j) = hyperbolicratio_error(testpars{j});
end
%err

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};

return;