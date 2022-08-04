%SS_TABULATE: organize spiking data from analyze_ss_responses into useful
%datasets.
%%% Last updated 11/04/15 on vivid

%% how to choose which 'spikes' variable to analyze (.spikes_raw vs. .spikes_cleaned)?

%%% Split responses into sf_fixed and tf_fixed, then average the identical
%%% stimulus sets. 
trials_tffixed = trialdata(trialdata.tf==ss_pars.tf_fixed,:);
trials_sffixed = trialdata(trialdata.sf==ss_pars.sf_fixed,:);
unique_tffixed = ss_pars.nAngles * ss_pars.n_diams * ss_pars.n_sfs; 
unique_sffixed = ss_pars.nAngles * ss_pars.n_diams * ss_pars.n_tfs;

% 'sps_unitmean' contains the spikes for each unit averaged over all
%     iterations
% 'sps_shankmean' contains spikes averaged over all sites and iterations
% 'sps_all' contains .... 
% the 'angleAvg' datasets contain data from the 'repAvg' datasets but
%     treating angles as identical; if repetitions per unique stim are
%     unequal, spikes are averaged over concatenated trials from all
%     angles, not [within angles then between angles]
if ss_pars.n_sfs > 0 
    repAvg_tffixed = dataset({NaN(unique_tffixed,5),'Angle','diam','sf','tf','sps_shankmean'},...
        {cell(unique_tffixed,1),'sps_all'},...
        {NaN(unique_tffixed,length(snipdata.unitnames)), 'sps_unitmean'},...
        'ObsNames',cellstr(num2str((1:unique_tffixed)')));
    angleAvg_tffixed = repAvg_tffixed(1:unique_tffixed/ss_pars.nAngles,:); angleAvg_tffixed.Angle=[];
end 
if ss_pars.n_tfs > 0 
    repAvg_sffixed = dataset({NaN(unique_sffixed,5),'Angle','diam','sf','tf','sps_shankmean'},...
        {cell(unique_sffixed,1),'sps_all'},...
        {NaN(unique_sffixed,length(snipdata.unitnames)), 'sps_unitmean'},...
        'ObsNames',cellstr(num2str((1:unique_sffixed)')));
    angleAvg_sffixed = repAvg_sffixed(1:unique_sffixed/ss_pars.nAngles,:); angleAvg_sffixed.Angle=[];
end

comboA = 0; comboB = 0; 
for Angle = 1:ss_pars.nAngles  
    for diam = 1:ss_pars.n_diams
        for sf = 1:ss_pars.n_sfs
            comboA = comboA + 1;
            match=trials_tffixed.sf==sf_vals(sf) &...
                  trials_tffixed.diam==diam_vals(diam) &...
                  trials_tffixed.Angle==angle_vals(Angle);
            repAvg_tffixed.sf(comboA) = sf_vals(sf);
            repAvg_tffixed.tf(comboA) = ss_pars.tf_fixed;
            repAvg_tffixed.diam(comboA) = diam_vals(diam);
            repAvg_tffixed.Angle(comboA) = angle_vals(Angle);
            repAvg_tffixed.sps_all(comboA) = {trials_tffixed.spikes(match,:)};
            repAvg_tffixed.sps_unitmean(comboA,:) = nanmean(trials_tffixed.spikes(match,:));
            repAvg_tffixed.sps_shankmean(comboA) = nanmean(repAvg_tffixed.sps_unitmean(comboA,:));
        end
        if ss_pars.n_tfs > 0 
            for tf = 1:ss_pars.n_tfs
                comboB = comboB + 1;
                match=trials_sffixed.tf==tf_vals(tf) &...
                      trials_sffixed.diam==diam_vals(diam) &...
                      trials_sffixed.Angle==angle_vals(Angle);
                repAvg_sffixed.sf(comboB) = ss_pars.sf_fixed;
                repAvg_sffixed.tf(comboB) = tf_vals(tf);
                repAvg_sffixed.diam(comboB) = diam_vals(diam);
                repAvg_sffixed.Angle(comboB) = angle_vals(Angle);
                repAvg_sffixed.sps_all(comboB) = {trials_sffixed.spikes(match,:)};
                repAvg_sffixed.sps_unitmean(comboB,:) = nanmean(trials_sffixed.spikes(match,:));
                repAvg_sffixed.sps_shankmean(comboB) = nanmean(repAvg_sffixed.sps_unitmean(comboB,:));
            end
        end
    end
end

% Get unit-averaged responses for each unit.
anglenames = strrep(cellfun(@(x)['Angle_',num2str(x)],num2cell(angle_vals),'UniformOutput',false),'.','_');
diamnames = strrep(cellfun(@(x)['Diam_',num2str(x)],num2cell(diam_vals),'UniformOutput',false),'.','_'); 
sf_names = strrep(cellfun(@(x)['SF_',num2str(x)],num2cell(sf_vals),'UniformOutput',false),'.','_');
tf_names = strrep(cellfun(@(x)['TF_',num2str(x)],num2cell(tf_vals),'UniformOutput',false),'.','_');
unitsByAngle = repmat({cell(length(snipdata.unitnames),1)},1,length(anglenames)+1);
unitnans = NaN(length(snipdata.unitnames),1);
allAngles_tffixed = dataset(unitnans, unitsByAngle{:},...
    'VarNames',[{'pref_angle' 'sps_pref_angle'} anglenames],'ObsNames', snipdata.unitnames);
if ss_pars.n_tfs > 0 
    allAngles_sffixed = allAngles_tffixed;
end

%% give option to declare that no preferred angle is found, use data from all angles

% Find preferred angle.
%%% preferred angle is computed using responses only the smallest diameter
%%%% need to add statistical test to determine whether or not cell actually
%%%% has preferered angle (?ttest of largest to next-largest or anova... or
%%%% von Mises, look at SD
angleTuningNans = NaN(ss_pars.nAngles,1);
nanCellForUnits = repmat({angleTuningNans},1,length(snipdata.unitnames));
meanSpsByAngle = dataset(nanCellForUnits{:},...
    'VarNames',{snipdata.unitnames{:}},'ObsNames',anglenames);
for angleInd = 1:ss_pars.nAngles
    for diamInd = 1:length(diamIndsForPrefAngle)
        match = trialdata.Angle == angle_vals(angleInd) &...
                trialdata.diam == diam_vals(diamInd);
        thisAngleSps = mean(trialdata.spikes(match,:)); % all spikes for all units to this angle at specified diams
        for unitInd = 1:length(snipdata.unitnames)
            meanSpsByAngle{angleInd,unitInd} = thisAngleSps(unitInd);
        end
    end
end

for chan = 1:length(snipdata.unitnames)
    for Angle = 1:ss_pars.nAngles
        if ss_pars.n_sfs > 0 
            allAngles_tffixed{snipdata.unitnames(chan),anglenames(Angle)} =...
                dataset({NaN(ss_pars.n_sfs,ss_pars.n_diams),diamnames{:}},'ObsNames',sf_names);
        end
        if ss_pars.n_tfs > 0 
            allAngles_sffixed{snipdata.unitnames(chan),anglenames(Angle)} =...
                dataset({NaN(ss_pars.n_tfs,ss_pars.n_diams),diamnames{:}},'ObsNames',tf_names); 
        end
        for diam = 1:ss_pars.n_diams
            if ss_pars.n_sfs > 0 
                for sf = 1:ss_pars.n_sfs
                    match = repAvg_tffixed.diam == diam_vals(diam) &...
                        repAvg_tffixed.sf == sf_vals(sf) &...
                        repAvg_tffixed.Angle == angle_vals(Angle);
                    allAngles_tffixed{snipdata.unitnames(chan),anglenames(Angle)}{sf_names(sf),diamnames(diam)} =...
                        repAvg_tffixed.sps_unitmean(match,chan);
                end
            end
            if ss_pars.n_tfs > 0 
                for tf = 1:ss_pars.n_tfs
                    match = repAvg_sffixed.diam == diam_vals(diam) &...
                        repAvg_sffixed.tf == tf_vals(tf) &...
                        repAvg_sffixed.Angle == angle_vals(Angle);
                    allAngles_sffixed{snipdata.unitnames(chan),anglenames(Angle)}{tf_names(tf),diamnames(diam)} =...
                        repAvg_sffixed.sps_unitmean(match,chan);
                end
            end
        end
    end
end
for unitInd = 1:length(snipdata.unitnames)
    [junk prefAngleInd] = max(double(meanSpsByAngle(:,unitInd))); % index of preferred angle 
    if ss_pars.n_sfs > 0 
        allAngles_tffixed.pref_angle(unitInd) = angle_vals(prefAngleInd);
        allAngles_tffixed.sps_pref_angle{unitInd} = allAngles_tffixed{unitInd,anglenames{prefAngleInd}}; % sps at preferred angle
    end
    if ss_pars.n_tfs > 0 
        allAngles_sffixed.pref_angle(unitInd) = angle_vals(prefAngleInd);
        allAngles_sffixed.sps_pref_angle{unitInd} = allAngles_sffixed{unitInd,anglenames{prefAngleInd}}; % sps at preferred angle
    end
end
    
%%%%% find best angle, fill in last 2 rows

% % % Check to see whether there were an equal number of repetitions for each
% % % unique stim parameter combination. (doesn't appear to be necessary)
% % tffixed_equalReps = 1; sffixed_equalReps = 1; % default
% % tffixed_nReps = cellfun(@(x)size(x,1),repAvg_tffixed.sps_all,'UniformOutput',false);
% % sffixed_nReps = cellfun(@(x)size(x,1),repAvg_sffixed.sps_all,'UniformOutput',false);
% % if length(unique(tffixed_nReps)) > 1;
% %     tffixed_equalReps = 0;
% % end
% % if length(unique(sffixed_nReps)) > 1;
% %     sffixed_equalReps = 0;
% % end

% Take average values over all angles. 
comboA = 0; comboB = 0;
for diam = 1:ss_pars.n_diams
    if ss_pars.n_sfs > 0 
        for sf = 1:ss_pars.n_sfs
            comboA = comboA+1;
            match = repAvg_tffixed.diam == diam_vals(diam) &...
                    repAvg_tffixed.sf == sf_vals(sf); % could limit this match search for speed
            angleAvg_tffixed.diam(comboA) = diam_vals(diam);
            angleAvg_tffixed.sf(comboA) = sf_vals(sf);
            angleAvg_tffixed.tf(comboA) = ss_pars.tf_fixed;

            %%% does the  commented-out step below do  anything sensible?
            %%%  also doesn't work with unequal repetitions
% %                 angleAvg_tffixed.sps_all(comboA) =... % average the sps_all values for all angles
% %                     {mean(reshape(cell2mat((repAvg_tffixed.sps_all(match))'),...
% %                     ss_pars.repetitions,length(snipdata.unitnames),ss_pars.nAngles),3)};
            angleAvg_tffixed.sps_all(comboA) = ... % avg sps_all for all identical-stim trials, ignoring angle
                {nanmean(cell2mat(repAvg_tffixed.sps_all(match)),1)}; % are some of this and the next 2 lines redundant?
            angleAvg_tffixed.sps_unitmean(comboA,:) = nanmean(angleAvg_tffixed.sps_all{comboA,:},1); % make sure mean is taken vertically if input is a single row
%             angleAvg_tffixed.sps_shankmean(comboA) = mean(angleAvg_tffixed.sps_unitmean(comboA,:),1);  %%% not yet computing shankmean... code may be wrong    
        end
    end
    if ss_pars.n_tfs > 0 
        for tf = 1:ss_pars.n_tfs
            comboB = comboB+1;
            match = repAvg_sffixed.diam == diam_vals(diam) &...
                    repAvg_sffixed.tf == tf_vals(tf); % could limit this match search for speed
            angleAvg_sffixed.diam(comboB) = diam_vals(diam);
            angleAvg_sffixed.sf(comboB) = ss_pars.sf_fixed;
            angleAvg_sffixed.tf(comboB) = tf_vals(tf);
            
            %%% does the  commented-out step below do  anything sensible?
            %%%  also doesn't work with unequal repetitions
% %                 angleAvg_sffixed.sps_all(comboA) =... % average the sps_all values for all angles
% %                     {mean(reshape(cell2mat((repAvg_sffixed.sps_all(match))'),...
% %                     ss_pars.repetitions,length(snipdata.unitnames),ss_pars.nAngles),3)};
            angleAvg_sffixed.sps_all(comboB) = ... % avg sps_all for all identical-stim trials, ignoring angle
                {nanmean(cell2mat(repAvg_sffixed.sps_all(match)),1)}; % are some of this and the next 2 lines redundant?
            angleAvg_sffixed.sps_unitmean(comboB,:) = nanmean(angleAvg_sffixed.sps_all{comboB,:},1); % make sure mean is taken vertically if input is a single row
%             angleAvg_sffixed.sps_shankmean(comboB) = mean(angleAvg_sffixed.sps_unitmean(comboB,:),1);  %%% not yet computing shankmean... code may be wrong    
        end
    end
end

% Construct site-specific size tuning curves averaged over angles and repetitions.
%%%% Each element of unitAvg_sffixed and unitAvg_tffixed contain data for
%%%% this unit; within the element for each unit, 'sps_all' contains number
%%%% of spikes responding to the diam and sf or tf listed for every
%%%% repetition, and 'sps_mean' lists the means of these trials.
unitcells = cell(length(unitnames),1);
if ss_pars.n_sfs > 0 
% % % %     unitAvg_tffixed = cell(length(snipdata.unitnames),3);
    unitAvg_tffixed = dataset(unitcells,unitcells,'ObsNames',unitnames,...
        'VarNames',{'sps_mean','sps_all'});
    sfdiamnan_ds = mat2dataset(NaN(ss_pars.n_diams,ss_pars.n_sfs),...
        'ObsNames',diamnames,'VarNames',sf_names); % for mean spikes
    sfdiamcells_ds = mat2dataset(cell(size(sfdiamnan_ds)),...
        'ObsNames',diamnames,'VarNames',sf_names); % for mean spikes
end
if ss_pars.n_tfs > 0 
% % % % %     unitAvg_sffixed = cell(length(snipdata.unitnames),3);
    unitAvg_sffixed = dataset(unitcells,unitcells,'ObsNames',unitnames,...
        'VarNames',{'sps_mean','sps_all'});
    tfdiamnan_ds = mat2dataset(NaN(ss_pars.n_diams,ss_pars.n_tfs),...
        'ObsNames',diamnames,'VarNames',tf_names); % for mean spikes
    tfdiamcells_ds = mat2dataset(cell(size(tfdiamnan_ds)),...
        'ObsNames',diamnames,'VarNames',tf_names); % for mean spikes
end

for unitInd = 1:length(snipdata.unitnames)
    thisUnit = unitnames(unitInd);
    if ss_pars.n_sfs > 0 
        unitAvg_tffixed{thisUnit,'sps_mean'} = sfdiamnan_ds; 
        unitAvg_tffixed{thisUnit,'sps_all'} = sfdiamcells_ds; 
        for diam = 1:ss_pars.n_diams % make sure sfreqs and tfreqs are taken in the correct order
            for sf = 1:length(sf_vals)
                match = angleAvg_tffixed.diam==diam_vals(diam) & angleAvg_tffixed.sf==sf_vals(sf);
                unitAvg_tffixed{thisUnit,'sps_mean'}{diamnames{diam},sf_names{sf}} = ...
                    angleAvg_tffixed.sps_unitmean(match,unitInd);
                match = repAvg_tffixed.diam==diam_vals(diam) & repAvg_tffixed.sf==sf_vals(sf);
                unitAvg_tffixed{thisUnit,'sps_all'}{diamnames{diam},sf_names{sf}} = ...
                    repAvg_tffixed.sps_all{match}(:,unitInd);
            end
        end
    end
    if ss_pars.n_tfs > 0 
        unitAvg_sffixed{thisUnit,'sps_mean'} = tfdiamnan_ds; 
        unitAvg_sffixed{thisUnit,'sps_all'} = tfdiamcells_ds; 
        for diam = 1:ss_pars.n_diams % make sure sfreqs and tfreqs are taken in the correct order
            for tf = 1:length(tf_vals)
                match = angleAvg_sffixed.diam==diam_vals(diam) & angleAvg_tffixed.tf==tf_vals(tf);
                unitAvg_sffixed{thisUnit,'sps_mean'}{diamnames{diam},tf_names{tf}} = ...
                    angleAvg_sffixed.sps_unitmean(match,unitInd);
                match = repAvg_sffixed.diam==diam_vals(diam) & repAvg_sffixed.tf==tf_vals(tf);
                unitAvg_sffixed{thisUnit,'sps_all'}{diamnames{diam},tf_names{tf}} = ...
                    repAvg_sffixed.sps_all{match}(:,unitInd);
            end
        end
    end    
end

% Construct size tuning curves averaged over angles, repetitions, and shanks. 
if ss_pars.n_sfs > 0 
    grandAvg_tffixed = dataset({NaN(ss_pars.n_sfs,ss_pars.n_diams),diamnames{:}},'ObsNames',sf_names);
end
if ss_pars.n_tfs > 0 
    grandAvg_sffixed = dataset({NaN(ss_pars.n_tfs,ss_pars.n_diams),diamnames{:}},'ObsNames',tf_names); 
end

for diam = 1:ss_pars.n_diams % make sure sfreqs and tfreqs are taken in the correct order
    if ss_pars.n_sfs > 0 
        for sf = 1:length(sf_vals)
            match = angleAvg_tffixed.diam==diam_vals(diam) & angleAvg_tffixed.sf==sf_vals(sf);
            grandAvg_tffixed.(diamnames{diam})(sf_names{sf}) = angleAvg_tffixed.sps_shankmean(match);
        end
    end
    if ss_pars.n_tfs > 0 
        for tf = 1:length(tf_vals) %%% not working? sometimes only see last row filled in
            match = angleAvg_sffixed.diam==diam_vals(diam) & angleAvg_sffixed.tf==tf_vals(tf);
            grandAvg_sffixed.(diamnames{diam})(tf_names{tf}) = angleAvg_sffixed.sps_shankmean(match);
        end
    end
end
