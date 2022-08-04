%NP_signif_test.m: This function takes in neurometric and psychometric data 
%	and calculates whether the slope and threshold are significantly different. 
%   Psychophysical data are assumed to behave like a binomial distribution, whereas
%   neurometric data are resampled with replacement. 95% confidence interval are calculated. 
%   Neurometric and psychometric data are significantly different if the 
%   confidence intervals do not overlap TU 08/13/02
function [P_alpha, P_beta, P_conf_alpha, P_conf_beta, N_conf_alpha, N_conf_beta] = Weibull_fit_signif_test(PC_psycho, PC_neuro, correct_trials, pref_dist, null_dist, X, N_trials)

% lets always make the inputs column vectors.  This way, either
% row or column vectors can be input and it will still work.
if ( (size(PC_psycho,2)>1) & (size(PC_psycho,1)==1) )	% a row vector, so transpose
   PC_psycho = PC_psycho';
end
if ( (size(PC_neuro,2)>1) & (size(PC_neuro,1)==1) )	% a row vector, so transpose
   PC_neuro = PC_neuro';
end
if ( (size(X,2)>1) & (size(X,1)==1) )	% a row vector, so transpose
   X = X';
end
if ( (size(N_trials,2)>1) & (size(N_trials,1)==1) )	% a row vector, so transpose
   N_trials = N_trials';
end

% loop through hundreds of times and do the following:
%   for psycho data
%	- randomly select a number from a binomial distribution and calculate percent correct
%	- fit the data with a weibull function
%   for neuro data 
%	- resample data with replacement and calculate ROC
%	- fit the data with a weibull function

MAX_ITER = 1000;
count_alpha = 0;
count_beta = 0;
for i=1:MAX_ITER
%    for j=1:length(X) %loop through each X value.
%        if ( N_trials(j) > 0)
%            %calculate percent correct
%            PC_temp1(j) = sum(binornd(1,PC_psycho(j),N_trials(j),1))/N_trials(j);  
%            %now, make data for Weibull fit
%            fit_data(j, 1) = X(j);
%            fit_data(j, 2) = PC_temp1(j);
%            fit_data(j, 3) = N_trials(j);
%        end
%    end
%    %now, fit a weibull function to PC_temp   
%    [P_bino_alpha(i) P_bino_beta(i)] = weibull_fit(fit_data);
%    
%    %do this for neuro data
%    for j=1:length(X) %loop through each X value.
%       if ( N_trials(j) > 0)
%            %calculate percent correct
%            PC_temp1(j) = sum(binornd(1,PC_neuro(j),N_trials(j),1))/N_trials(j);
%            %now, make data for Weibull fit
%            fit_data(j, 1) = X(j);
%            fit_data(j, 2) = PC_temp1(j);
%            fit_data(j, 3) = N_trials(j);
%        end
%    end
%    %now, fit a weibull function to PC_temp   
%    [N_bino_alpha(i) N_bino_beta(i)] = weibull_fit(fit_data);
    
    %resample data with replacement for psycho data
    for j=1:length(X) %loop through each X value.
        if ( N_trials(j) > 0)
            index = fix(rand(1,length(correct_trials{j}))*length(correct_trials{j})) + 1;
            correct_trials_temp = correct_trials{j}(index);
            %calculate percent correct
            PC_temp2(j) = sum(correct_trials_temp)/N_trials(j); 
            %now, make data for Weibull fit
            fit_data(j, 1) = X(j);
            fit_data(j, 2) = PC_temp2(j);
            fit_data(j, 3) = N_trials(j);
        end
    end
    [P_resamp_alpha(i) P_resamp_beta(i)] = weibull_fit(fit_data);
    
    %resample data with replacement for neuro data
    for j=1:length(X) %loop through each X value.
        if ( N_trials(j) > 0)
            pref_index = fix(rand(1,length(pref_dist{j}))*length(pref_dist{j})) + 1;
            pref_dist_temp = pref_dist{j}(pref_index);
            null_index = fix(rand(1,length(null_dist{j}))*length(null_dist{j})) + 1;
            null_dist_temp = null_dist{j}(null_index);
            %calculate percent correct
            PC_temp2(j) = rocN(pref_dist_temp, null_dist_temp, 100); 
            %now, make data for Weibull fit
            fit_data(j, 1) = X(j);
            fit_data(j, 2) = PC_temp2(j);
            fit_data(j, 3) = N_trials(j);
        end
    end                      
    [N_resamp_alpha(i) N_resamp_beta(i)] = weibull_fit(fit_data);
    
    diff_alpha(i) = P_resamp_alpha(i) - N_resamp_alpha(i);
    diff_beta(i) = P_resamp_beta(i) - N_resamp_beta(i);
end

%find confidence intervals
conf_interval = [2.5 97.5];
P_conf_alpha = PRCTILE(P_resamp_alpha, conf_interval);
N_conf_alpha = PRCTILE(N_resamp_alpha, conf_interval);
P_conf_beta = PRCTILE(P_resamp_beta, conf_interval);
N_conf_beta = PRCTILE(N_resamp_beta, conf_interval);

diff_conf_alpha = PRCTILE(diff_alpha, conf_interval);
diff_conf_beta = PRCTILE(diff_beta, conf_interval);

P_alpha = 0; P_beta = 0;
%if ((min(P_conf_alpha) > max(N_conf_alpha)) |  (min(N_conf_alpha) > max(P_conf_alpha)))
if ((min(diff_conf_alpha) > 0) | (0 > max(diff_conf_alpha)))
        P_alpha = 1;
end
%if ((min(P_conf_beta) > max(N_conf_beta)) |  (min(N_conf_beta) > max(P_conf_beta)))
if ((min(diff_conf_beta) > 0) | (0 > max(diff_conf_beta)))
        P_beta = 1;
end

%P_alpha_P = signrank(P_bino_alpha, P_resamp_alpha);
%P_alpha_N = signrank(N_bino_alpha, N_resamp_alpha);
%P_beta_P = signrank(P_bino_beta, P_resamp_beta);
%P_beta_N = signrank(N_bino_beta, N_resamp_beta);

return;