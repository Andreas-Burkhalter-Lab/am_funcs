%Weibull_fit_signif_test.m: This function takes in two vectors of percent correct data 
%	and calculates whether the slope and threshold are significantly different. Data are
%   assumed to behave like a binomial distribution.  TU 08/08/01
function [P_alpha, P_beta] = Weibull_fit_signif_test(PC1, PC2, X, N_trials)

% lets always make the inputs column vectors.  This way, either
% row or column vectors can be input and it will still work.
if ( (size(PC1,2)>1) & (size(PC1,1)==1) )	% a row vector, so transpose
   PC1 = PC1';
end
if ( (size(PC2,2)>1) & (size(PC2,1)==1) )	% a row vector, so transpose
   PC2 = PC2';
end
if ( (size(X,2)>1) & (size(X,1)==1) )	% a row vector, so transpose
   X = X';
end
if ( (size(N_trials,2)>1) & (size(N_trials,1)==1) )	% a row vector, so transpose
   N_trials = N_trials';
end

%fit a weibull function to PC1 and PC2
fit_data(:, 1) = X;
fit_data(:, 2) = PC1;
fit_data(:, 3) = N_trials;
[PC1_alpha PC1_beta] = weibull_fit(fit_data);
fit_data(:, 2) = PC2;
[PC2_alpha PC2_beta] = weibull_fit(fit_data);    

p_mean = (PC1+PC2)/2; %calculate the average of PC1 and PC2

% loop through hundreds of times and do the following:
%	- randomly select a number from a binomial distribution and calculate percent correct
%	- fit the data with a weibull function
%	- determine the proportion of times that the difference between the calculated parameters are greater than the difference between the original inputs

MAX_ITER = 2000;
count_alpha = 0;
count_beta = 0;
for i=1:MAX_ITER
    for j=1:length(PC1) %loop through each X value.
        if ( N_trials(j) > 0)
            %calculate percent correct
            PC_temp1(j) = sum(binornd(1,p_mean(j),N_trials(j),1))/N_trials(j);
            %now, make data for Weibull fit
            fit_data(j, 1) = X(j);
            fit_data(j, 2) = PC_temp1(j);
            fit_data(j, 3) = N_trials(j);
        end
    end
    %now, fit a weibull function to PC_temp   
    [temp1_alpha temp1_beta] = weibull_fit(fit_data);
    
    %Do it again
    for j=1:length(PC1) %loop through each X value.
        if ( N_trials(j) > 0)
            %calculate percent correct
            PC_temp2(j) = sum(binornd(1,p_mean(j),N_trials(j),1))/N_trials(j);
            %now, make data for Weibull fit
            fit_data(j, 1) = X(j);
            fit_data(j, 2) = PC_temp2(j);
            fit_data(j, 3) = N_trials(j);
        end
    end
    %now, fit a weibull function to PC_temp   
    [temp2_alpha temp2_beta] = weibull_fit(fit_data);
    
    if ( abs(temp1_alpha-temp2_alpha) > abs(PC1_alpha-PC2_alpha) )
        count_alpha = count_alpha + 1;
    end
    if ( abs(temp1_beta-temp2_beta) > abs(PC1_beta-PC2_beta) )
        count_beta = count_beta + 1;
    end
end

P_alpha = count_alpha/MAX_ITER;
P_beta = count_beta/MAX_ITER;

return;