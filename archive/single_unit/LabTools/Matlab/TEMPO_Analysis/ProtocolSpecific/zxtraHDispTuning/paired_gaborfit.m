function [pars] = paired_gaborfit(means,raw)
%means and raw are each cell arrays containing the data for each different group (ie, speed)

global Data1 Data2 RawData1 RawData2;

%allows calculation of error for both raw data and mean data
Data1 = means{1};
RawData1 = raw{1};
Data2 = means{2};
RawData2 = raw{2};

% first, generate some initial parameter guesses
[max_val_1	max_indx_1] = max(Data1(:,2));
[min_val_1	min_indx_1] = min(Data1(:,2));
[max_disp_1 max_disp_indx_1] = max(Data1(:,1));
[min_disp_1 min_disp_indx_1] = min(Data1(:,1));
[max_val_2	max_indx_2] = max(Data2(:,2));
[min_val_2	min_indx_2] = min(Data2(:,2));
[max_disp_2 max_disp_indx_2] = max(Data2(:,1));
[min_disp_2 min_disp_indx_2] = min(Data2(:,1));

N_values = length(Data1(:,1));

q(1) = mean(Data1(:,2));
q(2) = 0.5*(max_val_1 - min_val_1);
q(3) = 0.0;
q(4) = 0.5; % 
q(5) = 0.5; % 
q(6) = -1.5;
q(7) = mean(Data2(:,2));
q(8) = 0.5*(max_val_2 - min_val_2);

%get a good estimate of the starting frequency from the Fourier Transform
%for ease, we'll just take this from Data2
FFT_LENGTH = 128;
FT = abs(fft(Data2(:,2) - mean(Data2(:,2)), FFT_LENGTH));
FT2 = FT(1:length(FT)/2);
[max_FT_ampl max_FT_indx] = max(FT2);
%compute the disparity interval
dx = (max_disp_2 - min_disp_2)/(N_values-1);
freq = (max_FT_indx/length(FT2))/(2*dx);
q(5) = freq;

% do a quick search for a good starting phase
errors = [];
N=20;
for i=1:N
   q(6) = -pi + 2*pi*(i-1)/N;  
   errors(i) = paired_gaborerr(q);
end
[min_err min_err_indx] = min(errors);
q(6) = -pi + 2*pi*(min_err_indx-1)/N;

%set up the bounds arrays
A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
LB=[0;0;min_disp_1;0.1;0.3*freq;-2*pi;0;0];
UB=[1.5*max_val_1;1.35*(max_val_1 - min_val_1);max_disp_1;2*(max_disp_1 - min_disp_1);3.0*freq;2*pi;1.5*max_val_2;1.35*(max_val_2 - min_val_2)];

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');
%OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000);

N_reps = 20;
wiggle = 0.3;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    testpars{j} = fmincon('paired_gaborerr',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = paired_gaborerr(testpars{j});
end
%err

%now find best fit and return the parameters
[min_err min_err_indx] = min(err);
pars = testpars{min_err_indx};
return;