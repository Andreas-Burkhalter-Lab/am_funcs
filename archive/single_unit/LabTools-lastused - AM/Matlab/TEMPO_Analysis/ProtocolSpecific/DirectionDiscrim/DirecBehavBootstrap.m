% bootstrap to detect deviation of psychothreshold and neural threshold,
% expressed as what angle (+-difference) is thought to be significantly
% different

%-----------------------------------------------------------------------------------------------------------------------
function [avg_pct_correct, CI] = BehavBootstrap(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);   

TEMPO_Defs;		
Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%get the coherence
coherence = data.dots_params(DOTS_COHER, :, PATCH1);
unique_coherence = munique(coherence');

%now, select trials that fall between BegTrial and EndTrial
trials = 1:length(coherence);
%a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

trial_outcomes = logical(data.misc_params(OUTCOME, :) == CORRECT);

bootstp_num = 1000;
repetition = floor(length(trial_outcomes)/length(unique_coherence));
for b = 1 : bootstp_num
    % bootstrap dataset first
    for k = 1 : length(unique_coherence)
        select_boot = logical( coherence == unique_coherence(k));
        behav_select = trial_outcomes(select_boot); % behavior data
        for j = 1 : repetition
            behav_select = behav_select( randperm(length(behav_select)) ); % behavior data
            behav_bootstrap(j) = behav_select(1);
        end
        psycho_correct_boot(b,k) = sum(behav_bootstrap) / length(behav_bootstrap);
    end
end

%calculate confidence interval
avg_pct_correct = mean(mean(psycho_correct_boot))
s = avg_pct_correct * (1-avg_pct_correct);
t = tinv(0.05/2, length(trial_outcomes));
CI = s*t/sqrt(length(trial_outcomes)-1)
CI_lb = avg_pct_correct + CI
CI_ub = avg_pct_correct - CI

%histogram of performances in the bootstraps
figure
m=sort(mean(psycho_correct_boot,2));
hist(m,ceil(100*range(m)))
[h,x] = hist(m,ceil(100*range(m)));
bootpct_lb = m(round(bootstp_num*0.025))
bootpct_ub = m(round(bootstp_num*0.975))
xlabel('% Correct'); ylabel('# Bootstraps');

    
%     % calculate ROC
%     fit_data_neuro_boot = [];
%     for k = 1 : length(unique_stim_type)
%         if i < (1+length(unique_heading))/2
%             Neuro_correct_boot{k}(i) =  rocN( resp_heading_boot{k}(i,:),resp_heading_boot{k}((1+length(unique_heading))/2,:),100 ); % compare to the 0 heading condition, which is straght ahead
%             Neuro_correct_boot{k}(i) = 1 - Neuro_correct_boot{k}(i); % turn proportion correct into rightward choice
%         else
%             Neuro_correct_boot{k}(i) =  rocN( resp_heading_boot{k}((1+length(unique_heading))/2,:), resp_heading_boot{k}(i+1,:),100 ); % compare to the 0 heading condition, which is straght ahead
%         end
%         if  resp_mat{k}(1) < resp_mat{k}(end)  
%             % we don't know whether left heading is for sure smaller than right heading,thus ROC curve might be flipped above or below the unity line
%             % here we asume if left response larger than right response then asign the left to be preferred direction
%             Neuro_correct_boot{k}(i) = 1 - Neuro_correct_boot{k}(i);            
%         end         
%     end
%     % fit gaussian
%     for k = 1 : length(unique_stim_type)
%         neu_heading = unique_heading(unique_heading~=0);
%         beta = [0, 1.0];
%         [betafit_ne_boot{k},resids_ne_boot{k},J_ne_boot{k}] = nlinfit(neu_heading, Neuro_correct_boot{k}(:), 'cum_gaussfit', beta);  % fit data with least square
%         neu_thresh_boot(k,b) = betafit_ne_boot{k}(2);
%         [betafit_psy_boot{k},resids_psy_boot{k},J_psy_boot{k}] = nlinfit(unique_heading, psycho_correct_boot(k,:), 'cum_gaussfit', beta);  % fit data with least square
%         psy_thresh_boot(k,b) = betafit_psy_boot{k}(2);  
%     end
% end
% % test confidence field
% bin_num = 100; % temporally set 100 bins
% for k = 1 : length(unique_stim_type)
%     hist_ne_boot(k,:) = hist( neu_thresh_boot(k,:), bin_num );  % for bootstrap
%     bin_ne_sum = 0;
%     n_ne = 0;
%     while ( bin_ne_sum < 0.05*sum( hist_ne_boot(k,:)) )   % define confidential value to be 0.05
%           n_ne = n_ne+1;
%           bin_ne_sum = bin_ne_sum + hist_ne_boot(k, n_ne) + hist_ne_boot(k, bin_num-n_ne+1);      
%           neu_boot(k) = betafit_ne{k}(2) - min(neu_thresh_boot(k,:)) - n_ne * ( (max(neu_thresh_boot(k,:))-min(neu_thresh_boot(k,:))) / bin_num) ;    % calculate what value is thought to be significant different
%     end 
%     hist_psy_boot(k,:) = hist( psy_thresh_boot(k,:), bin_num );  % psycho data
%     bin_psy_sum = 0;
%     n_psy = 0;
%     while ( bin_psy_sum < 0.05*sum( hist_psy_boot(k,:)) )   % define confidential value to be 0.05
%           n_psy = n_psy + 1;
%           bin_psy_sum = bin_psy_sum + hist_psy_boot(k, n_psy) + hist_psy_boot(k,bin_num-n_psy+1);      
%           psy_boot(k) = betafit{k}(2) - min(psy_thresh_boot(k,:)) - n_psy * ( (max(psy_thresh_boot(k,:))-min(psy_thresh_boot(k,:))) / bin_num);   
%     end 
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%