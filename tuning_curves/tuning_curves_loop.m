%%%% get tuning curves for each stimulus parameter  tested
% updated 2020/4/20 on thermaltake


for thispar = pars.stimpars_to_analyze
    do_analysis = 0; % default
    thispar = thispar{:};
    rounded_diams = 1e-6 * round(1e6*stimpar_sets.diam); % sometimes diam vals are shifted slightly; round back to original vals
    initialize_stim_parameter(); 
    if do_analysis
        % for each roi, make table of F response for each trial sorted by stim parameter and trial 
        tuning = table(NaN(nrois,2),false(nrois,1),NaN(nrois,1),false(nrois,1),NaN(nrois,1),NaN(nrois,1),     cell(nrois,3),          cell(nrois,3),     cell(nrois,1),   cell(nrois,1),  NaN(nrois,2),       cell(nrois,1),...
            'VariableNames',{'centeryx','sgnf','anovap',    'sgnf_rspv',    'anovap_rspv',    'best_tested_val','timecourse_mean_best','timecourse_sem_best','trials','roi_image_prereg','centeryx_pre_rotate','roi_image_pre_rotate'});
        wbar = waitbar(0,['Getting ' thispar ' tuning curves...']);
        for iroi = 1:nrois 
            if isvalid(wbar)
                waitbar(iroi/nrois,wbar)
            end
            
            %%%% record the xy location and center of an roi into the tuning table and rotate if appropriate               
            record_roi_location();
            
            % get responses, test for significance; 'timecourse' = trial-averaged dff timecourse of response to each stim paramater value
            %   'timecourse' variables are split into prestim, during-stim, and poststim epochs
            tuning.trials{iroi} = table(stimparvals,NaN(nparvals,1),cell(nparvals,3),cell(nparvals,3),NaN(nparvals,stimpars.repetitions),NaN(nparvals,stimpars.repetitions),...
                'VariableNames',{       thispar,   'resp_mean',     'timecourse_mean', 'timecourse_sem',    'resp',                             'prestim_base'});
            for indpar = 1:nparvals % get dff values during stim and during prestim window
                match = find(stimpar_sets{:,par_sets_varname}==stimparvals(indpar) & otherparsfixed);
                sqrt_match = sqrt(length(match)); % for standard error of mean
                tuning.trials{iroi}.resp(indpar,1:length(match)) = [stimpar_sets.dff_during_stim(match,iroi)]';
                tuning.trials{iroi}.resp_mean(indpar) = nanmean(tuning.trials{iroi}.resp(indpar,:));
                tuning.trials{iroi}.prestim_base(indpar,1:length(match)) = [stimpar_sets.dff_prestim(match,iroi)]';
                prestim_timecourse_all_matching_trials = cell2mat(stimpar_sets.timecourse_prestim(match,iroi)')'; 
                stim_timecourse_all_matching_trials = cell2mat(stimpar_sets.timecourse_stim(match,iroi)')';
                poststim_timecourse_all_matching_trials = cell2mat(stimpar_sets.timecourse_poststim(match,iroi)')';
                % save mean and standard error timecourses
                tuning.trials{iroi}.timecourse_mean{indpar,1} = nanmean(prestim_timecourse_all_matching_trials);
                tuning.trials{iroi}.timecourse_mean{indpar,2} = nanmean(stim_timecourse_all_matching_trials);
                tuning.trials{iroi}.timecourse_mean{indpar,3} = nanmean(poststim_timecourse_all_matching_trials);
                tuning.trials{iroi}.timecourse_sem{indpar,1} = nanstd(prestim_timecourse_all_matching_trials) ./ sqrt_match;
                tuning.trials{iroi}.timecourse_sem{indpar,2} = nanstd(stim_timecourse_all_matching_trials)./ sqrt_match;
                tuning.trials{iroi}.timecourse_sem{indpar,3} = nanstd(poststim_timecourse_all_matching_trials) ./ sqrt_match;
            end
            tuning.anovap(iroi) = anova1(tuning.trials{iroi}.resp',[],'off'); % transpose resp so that responses are grouped by stim parameter value
            tuning.sgnf(iroi) = tuning.anovap(iroi) < 0.05;
            
            %%% commented out because this value gets overwritten by other method for computing responsivity
% % % % % % % % % % %             %%%% test for responsivity by performing anova on different
% % % % % % % % % % %             %%%% param vals plus baseline activity
% % % % % % % % % % %             resp_for_anova = tuning.trials{iroi}.resp'; % transpose so that resps will be extracted in order of param val, not trial number
% % % % % % % % % % %             resp_for_anova = resp_for_anova(:);
% % % % % % % % % % %             base_for_anova = tuning.trials{iroi}.prestim_base'; % transpose so that resps will be extracted in order of param val, not trial number
% % % % % % % % % % %             base_for_anova = base_for_anova(:);
% % % % % % % % % % %             dff_for_anova = [resp_for_anova; base_for_anova];
% % % % % % % % % % %             stim_names_for_anova = cellstr(num2str(kron(tuning.trials{iroi}{:,thispar},ones(size(tuning.trials{iroi}.resp,2),1)))); % dummy vars for different stim params
% % % % % % % % % % %             base_names_for_anova = repmat({'dff_prestim'},length(base_for_anova),1); % dummy vars for base
% % % % % % % % % % %             combonames_for_anova = [stim_names_for_anova; base_names_for_anova];
% % % % % % % % % % %             tuning.anovap_rspv(iroi) = anova1(dff_for_anova, combonames_for_anova, 'off'); % get responsiveness anova p value
% % % % % % % % % % %             tuning.sgnf_rspv(iroi) = tuning.anovap_rspv(iroi) < 0.05;
            
            % fit functions to tuning curves 
            %%%% hwhm_spacing = intervals to check when looking for preferred stim param value
            %%%% hwhm xmax = max value that we bother to compute as preferred... larger values take longer
            if pars.fit_tuning_functions
                meanforfit = [stimparvals, tuning.trials{iroi}.resp_mean];
                respforfit = [repmat(stimparvals,stimpars.repetitions,1), tuning.trials{iroi}.resp(:)]; % first col = stimparval, second = response val
                respforfit( isnan(respforfit(:,2)),: ) = []; % delete nan trials before fitting
                switch thispar
                    case 'sf'
                        fitparams = SF_loggaussfit(meanforfit,respforfit); %%% fitting function from Gao et al. 2011
                        fitfunc = @(q,stimval)q(2).*exp(-1./2./q(4).^2.*(log((stimval+q(5))./(q(3)+q(5)))).^2)+q(1); % lognormal from Gao et al. 2010
                        tuning.sf_pref(iroi) = fitparams(3); 
                        tuning.sf_fit_base_F(iroi) = fitparams(1);
                        tuning.sf_fit_amp(iroi) = fitparams(2); % amplitude
                        tuning.sf_fit_sigma(iroi) = fitparams(4); % standard deviation
                        tuning.sf_fit_log_offset(iroi) = fitparams(5);
                        tuning.sf_fitparams{iroi} = fitparams;
                        tuning.Properties.UserData = 'sf in cycles per degree';
                        prefval = tuning.sf_pref(iroi);
                        hwhm_xspacing = 0.0005;
                        hwhm_xmax = 3.5; 
                    case 'tf'
                        fitparams = TF_loggaussfit(meanforfit,respforfit); %%% fitting function from Gao et al. 2011
                        fitfunc = @(q,stimval)q(2).*exp(-1./2./q(4).^2.*(log((stimval+q(5))./(q(3)+q(5)))).^2)+q(1); % lognormal from Gao et al. 2010   
                        tuning.tf_pref(iroi) = fitparams(3);
                        tuning.tf_fit_base_F(iroi) = fitparams(1);
                        tuning.tf_fit_amp(iroi) = fitparams(2); % amplitude
                        tuning.tf_fit_sigma(iroi) = fitparams(4); % standard deviation
                        tuning.tf_fit_log_offset(iroi) = fitparams(5);
                        tuning.tf_fitparams{iroi} = fitparams;
                        tuning.Properties.UserData = 'tf in hertz';
                        prefval = tuning.tf_pref(iroi);
                        hwhm_xspacing = 0.01; 
                        hwhm_xmax = 30;  
                    case 'orient'
                        fitparams = von_mises_fit(meanforfit,respforfit); %%% fitting function from Gao et al. 2011
                        fitfunc = @(q,stimval)q(1).* exp( q(2).* (cos(stimval.*pi/180-q(3)*pi/180)-1)) + q(4).* exp( q(5).* (cos(stimval.*pi/180-q(6)*pi/180)-1)) + q(7);  % von mises from Gao et al. 2010                    
                        if fitfunc(fitparams,fitparams(3)) < fitfunc(fitparams,fitparams(6)) % if the second pref orient generates higher response than first pref orient, switch their order
                            fitparams = [fitparams(4:6); fitparams(1:3); fitparams(7)]; % put most preferred orient first
                        end
                        tuning.orient_pref1(iroi) = mod(fitparams(3), 360); % make sure pref orient is between 0 and 360
                        tuning.orient_amp1(iroi) = fitparams(1);
                        tuning.orient_width1(iroi) = fitparams(2);
                        tuning.orient_pref2(iroi) = mod(fitparams(6), 360); % make sure pref orient is between 0 and 360
                        tuning.orient_amp2(iroi) = fitparams(4);
                        tuning.orient_width2(iroi) = fitparams(5);
                        tuning.orient_fit_base_F(iroi) = fitparams(7);
                        tuning.orient_fitparams{iroi} = fitparams;
                        orient_pref1_ortho = tuning.orient_pref1(iroi) + 90; % least preferred orient; see selectivity index, OSI from Hofer et al. 2011
                        response_best = fitfunc(fitparams,tuning.orient_pref1(iroi)); % response_best from Hofer et al. 2011
                        response_ortho = fitfunc(fitparams,orient_pref1_ortho); % response_ortho from Hofer et al. 2011
                        tuning.orient_si(iroi) = [response_best - response_ortho] / [response_best + response_ortho]; % OSI from Hofer et al. 2011
                        tuning.Properties.UserData = 'orient in degrees';
                        prefval = tuning.orient_pref1(iroi);
                        hwhm_xspacing = 0.1;
                        hwhm_xmax = tuning.orient_pref1(iroi) + 200;
                    case 'diam'
                        fitparams = diff_erf_fit(meanforfit,respforfit);  %%% fitting function from Gao et al. 2011
                        fitfunc = @(q,stimval)q(1)*erf(stimval/q(2))-q(3)*erf(stimval/(q(2)+q(4)))+q(5);     %%% difference-of-error function plus baseline
                        tuning.diam_excit_size(iroi) = fitparams(2); % size of the excitatory center
                        tuning.diam_excit_amp(iroi) = fitparams(1);
                        tuning.diam_inhib_size(iroi) = fitparams(4);% size of the inhibitory surround
                        tuning.diam_inhib_amp(iroi) = fitparams(3);
                        tuning.diam_fit_base_F(iroi) = fitparams(5);
                        tuning.Properties.UserData = 'stimulus diameter in degrees';
                        prefval = tuning.diam_excit_size(iroi);
                        hwhm_xspacing = 0.05; 
                        hwhm_xmax = 100; 
                end
                % search from response peak to hwhm_xmax for the right-side hwhm  
                if prefval+hwhm_xspacing >= hwhm_xmax
                    tuning{iroi,[thispar, '_hwhm']} = NaN; % tuning curve is so broad, don't bother computing hwhm
                else
                    tuning{iroi,[thispar, '_hwhm']} = hwhm(fitfunc, fitparams, [prefval  hwhm_xmax], hwhm_xspacing, prefval, tuning{iroi,[thispar,'_fit_base_F']}); 
                end
                % use prefval to get nearest tested parval, save timecoures mean and sem to tuning table for that val
                [~,bestind] = min(abs(stimparvals-prefval)); % could compute difference on log scale
                tuning.best_tested_val(iroi) = stimparvals(bestind); %% closest tested stimparval to the fitted preferred value
                tuning.timecourse_mean_best(iroi,:) = tuning.trials{iroi}.timecourse_mean(bestind,:); % timecourse of mean F of best tested stimparval
                tuning.timecourse_sem_best(iroi,:) = tuning.trials{iroi}.timecourse_sem(bestind,:);  % timecourse of standard error of F of best tested stimparval
            end    
        end   
        if isvalid(wbar)
            close(wbar)
        end    

         %%%%% consolidate tuning_tables
        tuning.Properties.VariableNames{strcmp(tuning.Properties.VariableNames,'trials')} = [thispar '_trials']; % rename var
        tuning.Properties.VariableNames{strcmp(tuning.Properties.VariableNames,'anovap')} = [thispar '_anovap']; % rename var
        tuning.Properties.VariableNames{strcmp(tuning.Properties.VariableNames,'anovap_rspv')} = [thispar '_anovap_rsvp']; % rename var
        tuning.Properties.VariableNames{strcmp(tuning.Properties.VariableNames,'sgnf_rspv')} = [thispar  '_sgnf_rspv']; % rename var
        tuning.Properties.VariableNames{strcmp(tuning.Properties.VariableNames,'sgnf')} = [thispar '_sgnf']; % rename var
% % % % % % % % % % % % % % % % % % % % %         tuning.Properties.VariableNames{strcmp(tuning.Properties.VariableNames,'si')} = [thispar '_si']; % rename var... selectivity index, OSI from Hofer et al. 2011
        tuning.Properties.VariableNames{strcmp(tuning.Properties.VariableNames,'best_tested_val')} = [thispar '_best_tested_val']; % rename var
        tuning.Properties.VariableNames{strcmp(tuning.Properties.VariableNames,'timecourse_mean_best')} = [thispar '_timecourse_mean_best']; % rename var
        tuning.Properties.VariableNames{strcmp(tuning.Properties.VariableNames,'timecourse_sem_best')} = [thispar '_timecourse_sem_best']; % rename var
        if isfield(res,'tuning') % if another stimpar has already been processed
            delvar = [];
            for i = 1:size(tuning,2)
                if iscell(tuning{1,i}) && ~ischar(tuning{1,i}) && any(strcmp(tuning.Properties.VariableNames{i},res.tuning.Properties.VariableNames))
                    delvar = [delvar i];  % non-char cell can't be a key var, so delete
                end
            end
            tuning(:,delvar) = [];
            tuning.centeryx_pre_rotate = []; %%%% this var might be NaNs
            res.tuning = join(res.tuning,tuning); % join previously processed stimpar with current stimpar in tuning table; keyvars include 'centeryx','rspv','rspv_pval'
        else
            res.tuning = tuning;
        end

        clear tuning stimparvals nparvals par_sets_varname otherparsfixed match
    end
end

