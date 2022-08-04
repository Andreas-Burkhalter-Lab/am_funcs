%% this version doesn't allow for as many variations of exclusions as later version (for example, exclude based on RF tuning)

%%% get pairwise correlations between cells
%%%%%%%% [cortable, dff_downsampled] = pairwise_cor(tuningdat, starttime_voltage_recording, endtime_voltage_recording, show_plots, pars)
%
% tuningdat = output from get_tuning_curves; must contain tuning table containing quantile variable
% starttime_voltage_recording = starting value in scope_events.onset from which to start analyze timepoints; if empty, will use earliest
% endtime_voltage_recording = ending value in scope_events.onset from which to stop analyze timepoints; if empty, will use latest
%%%%% updated 19/5/3 on thermaltake

function [cortable, dff_downsampled] = pairwise_cor(tuningdat, starttime_voltage_recording, endtime_voltage_recording, show_plots, pars)

only_use_tuned_or_responsive_rois = 0; % only analyze cells that are tuned for sf, tf, or orientation or are responsive to stimuli
compute_noise_correlations = 1; 
    pars.loc_range_mps = field_default(pars,'loc_range_mps',[-inf inf]); %%% for noise correlations, only include trials with mean locomotion in meters/s within this range
stimparams = {'sf','tf','orient',}; %%% parameters to analyze for noise correlations - don't include rf
show_plots = vardefault('show_plots',0); 
% downsample dff (by averaging) before computing correlations to see how sample rate affects correlations
pars.downsample_dff = field_default(pars,'downsample_dff',0); 
    pars.downsample_factor = field_default(pars,'downsample_factor',1);

nrois = height(tuningdat.tuning);
ncombos = ((nrois-1)^2 + nrois-1) / 2; % from (n^2 + n)/2 replacing n with n-1... sum of consecutive integers
% patchdata = load(tuningdat.patchdata_file); patchdata = patchdata.patchdata;
patchdata = tuningdat.patchdata;
umPerPixels = umPerPix(patchdata.zoom,patchdata.scope);



%%%% make cortable:
%      roi1 and roi2 = row from tuningdat.tuning from which the two rois were taken
%      corrcoef = Pearson correlation coefficient between the dff traces of the rois derived from dff matrix
%      umdist = xy micron distance between the rois measured in this plane
%      quantdist = difference in quantiles between the rois, as found by findpatches.m
cortable = table(NaN(ncombos,1),NaN(ncombos,1),NaN(ncombos,1),NaN(ncombos,1),NaN(ncombos,1),NaN(ncombos,1),NaN(ncombos,1),NaN(ncombos,2),NaN(ncombos,2),...
    'VariableNames',{'roi1',     'roi2',     'quantile1',       'quantile2',       'corcoef_rest',     'umdist',    'quantdist',    'centeryx1', 'centeryx2'});

% find rows of scope_events to analyze, get correlation coefficients
if ~exist('starttime_voltage_recording','var') || isempty(starttime_voltage_recording) || isnan(starttime_voltage_recording)
    starttime_voltage_recording = min(tuningdat.scope_events.onset); % if not specified, start at beginning of scope_events
end
if ~exist('endtime_voltage_recording','var') || isempty(endtime_voltage_recording) || isnan(endtime_voltage_recording)
    endtime_voltage_recording = max(tuningdat.scope_events.onset); % if not specified, stop at end of scope_events
end
event_rows_to_analyze = find(tuningdat.scope_events.onset >= starttime_voltage_recording & tuningdat.scope_events.onset <= endtime_voltage_recording);

% get trials which exceed the locomotion threshold
locomotion_keep_trials = struct;
for iparam = 1:length(stimparams)
    thispar = stimparams{iparam};
    locomotion_keep_trials.(thispar) = [tuningdat.tuningpars.locm_trials.(thispar).locm_forw_mps > pars.loc_range_mps(1) & ...
                                        tuningdat.tuningpars.locm_trials.(thispar).locm_forw_mps < pars.loc_range_mps(2)]';
end

% get corr coefs, downsampling if specified
raw_dff = tuningdat.scope_events.dff(event_rows_to_analyze,:); % dff before downsampling
if pars.downsample_dff
    remaining_timepoints = mod(size(raw_dff,1), pars.downsample_factor); % timepoints not dividing evenly into downsample_factor
    raw_dff = raw_dff(1:end-remaining_timepoints,:); % elimate timepoints at end to make averaging easier
    % rearrange so we can average/downsample
    rearranged_dff_raw = permute(reshape(raw_dff,pars.downsample_factor,round(size(raw_dff,1)/pars.downsample_factor),[]), [2 3 1]); 
    dff_mat = mean(rearranged_dff_raw,3); % average along z to downsample/average every downsample_factor timepoints
    dff_downsampled = dff_mat;
else % don't downsample
    dff_mat = raw_dff;
end
cor_coefs = corrcoef(dff_mat); % get coefficients
if any(ismember(tuningdat.tuning.Properties.VariableNames,'rspv_pval'))
    tuned_or_responsive_rows = [tuningdat.tuning.sf_sgnf | tuningdat.tuning.tf_sgnf | tuningdat.tuning.orient_sgnf | tuningdat.tuning.rspv_pval < .05]; 
elseif ~any(ismember(tuningdat.tuning.Properties.VariableNames,'rspv_pval'))
    tuned_or_responsive_rows = [tuningdat.tuning.sf_sgnf | tuningdat.tuning.tf_sgnf | tuningdat.tuning.orient_sgnf]; 
end
% tabulate roi pairs
icombo = 0; % counter to go down cortable rows
wbar = waitbar(icombo/ncombos,'Computing dff correlations...');
for roi1 = 1:nrois
    for roi2 = roi1+1:nrois % pair roi1 with every roi it hasn't yet been paired with
        if ~only_use_tuned_or_responsive_rois || [tuned_or_responsive_rows(roi1) && tuned_or_responsive_rows(roi2)]
            icombo = icombo+1;
            cortable.roi1(icombo) = roi1;
            cortable.roi2(icombo) = roi2;
            cortable.centeryx1(icombo,:) = tuningdat.tuning.centeryx_reg(roi1,:);
            cortable.centeryx2(icombo,:) = tuningdat.tuning.centeryx_reg(roi2,:);
            cortable.quantile1(icombo) = tuningdat.tuning.quantile(roi1);
            cortable.quantile2(icombo) = tuningdat.tuning.quantile(roi2);
            cortable.quantdist(icombo) = abs(cortable.quantile1(icombo) - cortable.quantile2(icombo)); % quantdist = difference in quantiles between the rois, as found by findpatches.m
            distpix = sqrt( [cortable.centeryx1(icombo,1)-cortable.centeryx2(icombo,1)]^2 + [cortable.centeryx1(icombo,2)-cortable.centeryx2(icombo,2)]^2 ); % euclidean pixel distance
            cortable.umdist(icombo) = distpix * umPerPixels; %      umdist = xy micron distance between the rois measured in this plane
            cortable.corcoef_rest(icombo) = cor_coefs(roi1,roi2); % stim-off correlation
         
            %% noise (stim-on) correlations.... could limit window from which stim-response values are taken to shorter 'transient response' period
            if compute_noise_correlations
                noisecors = []; % initialize list of noise cors for all stimpars and stimpar values
                for iparam = 1:length(stimparams) %%% get noise correlations for each stim param tested, then combine them
                    thispar = stimparams{iparam};
                    n_par_values = tuningdat.stimpars.(['n_', thispar, 's']);                    
                    if n_par_values > 0 %%% if this stimpar was tested for
                        % next 2 lines might be unnecessarily slow - call directly from struct rather than copying
                        resps_1 = tuningdat.tuning{roi1,[thispar,'_trials']}{:}.resp'; % responses of this roi to all tests of this stimpar, with trials running down rows, param val across columns
                        resps_2 = tuningdat.tuning{roi2,[thispar,'_trials']}{:}.resp'; % responses of this roi to all tests of this stimpar, with trials running down rows, param val across columns
                        % eliminate trials which don't exceed locomotion threshold
                        resps_1(~locomotion_keep_trials.(thispar)) = NaN;
                        resps_2(~locomotion_keep_trials.(thispar)) = NaN;
                        colcors = NaN(1,n_par_values); % cor for each param val
                        for i_par_val = 1:n_par_values % for each value of this stimpar, get the pairwise noise correlations
                            thiscolcor = corrcoef(resps_1(:,i_par_val), resps_2(:,i_par_val), 'Rows', 'complete'); % cor across this param val (down this column within resps)
                            colcors(i_par_val) = thiscolcor(2,1);
                        end
                        noisecors = [noisecors, colcors]; % add correlations from this stimpar to the complete list
                        cortable.corcoef_noise(icombo) = nanmean(noisecors); % noise correlations averaged over every stim value from every stimpar
                    end
                end
            end
            %%
        end
    end
    try wbar = waitbar(icombo/ncombos,wbar); end
end
try close(wbar); end
deleterows = isnan(cortable.roi1);
cortable(deleterows,:) = [];

%% plotting
if show_plots
    pairwise_cor_plotting(cortable);
end
