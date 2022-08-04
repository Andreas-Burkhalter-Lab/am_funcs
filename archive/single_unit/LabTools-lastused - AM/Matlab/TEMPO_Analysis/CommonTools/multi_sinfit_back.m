function [pars] = multi_sinfit(means,raw, q)
global Data RawData Data1 RawData1 Data2 RawData2 Data3 RawData3;

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*' 'c*'};
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--' 'c--'};

Data = cell(size(means));
RawData = cell(size(raw));

for i=1:length(means)
   Data{i} = means{i};
   RawData{i} = raw{i};
end

for i=1:length(Data)
   [max_val_array(i)	max_indx_array(i)] = max(Data{i}(:,2));
   [min_val_array(i)	min_indx_array(i)] = min(Data{i}(:,2));
   [max_disp(i) max_disp_indx(i)] = max(Data{i}(:,1));
   [min_disp(i) min_disp_indx(i)] = min(Data{i}(:,1));
   N_values(i) = length(Data{i}(:,1));
end

min_val = min(min_val_array);
max_val = min(max_val_array);


%initialize parameters--------------------------------------------------
%amplitude for data1
%q(1) = (max_val_array(1)-min_val_array(1))/2;

%frequency
%q(2) = 1;

%baseline for data1
%q(4) = mean(Data{1}(:,2));

%exponential
%q(5) = 1;

%index = 6;
%for i=2:length(Data)
%   q(index) = (max_val_array(i)-min_val_array(i))/2; %amp
%   q(index+1) = mean(Data{i}(:,2)); %baseline
%   q(index+2) = 1; %exp
%   index = index + 3;
%end

% do a quick search for a good starting phase
%errors = [];
%N=20;
%for i=1:N
%   q(3) = -pi + 2*pi*(i-1)/N;  
%   errors(i) = multi_sinerr(q);
%end
%[min_err min_indx] = min(errors);
%q(3) = -pi + 2*pi*(min_indx-1)/N;

% do a quick search for a good starting exp
%errors = [];
%N=40;
%rand_list=rand(N,length(Data));
%rand_list(2:2:length(rand_list),:) = rand_list(2:2:length(rand_list), :);

%for i=1:N
%   q(5) = rand_list(i,1);
%   index = 6;
%   for j=2:length(Data)
%      q(index+2) = rand_list(i,j);
%      index = index + 3;
%   end   
%   errors(i) = multi_sinerr(q);
%end
%[min_err min_indx] = min(errors);
%q(5) = rand_list(min_indx,1);
%index = 6;
%for j=2:length(Data)
%   q(index+2) = rand_list(min_indx,j);
%   index = index + 3;
%end   
%initialize parameters--------------------------------------------------

A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=[ 0; .6; -2*pi; min_val_array(1); -2];
UB=[ (max_val_array(1)-min_val_array(1))*2; 1.4; 2*pi; max_val_array(1); 2];

for i=2:length(Data)
    LB = [LB; 0; min_val_array(i); -2];
    UB = [UB; (max_val_array(i)-min_val_array(i))*2; max_val_array(i); 2];
end

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'on', 'MaxFunEvals', 2500);

N_reps = 40;
wiggle = 0.2;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    pause(.5);
    testpars{j} = fmincon('multi_sinerr',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = multi_sinerr(testpars{j});
    j
 end

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};
return;