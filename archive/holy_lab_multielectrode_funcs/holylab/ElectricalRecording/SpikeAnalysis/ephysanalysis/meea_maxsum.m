function analysis_out = meea_maxsum(ephysin, analysis_in)
% meea_maxsum chooses finds tmax and rmax for multielectrode ephys analysis
%
% syntax: analysis_out = meea_maxsum_tstart(ephysin, analysis_in)
%
% inputs:  ephysin: multi electrode ephys structure
%          analysis_in: analysis structure with optional subfields:
%            .trange:             {
%            .resample_rate:     {{ if not present, will be set to defaults of 'meea_resample'
%            .resampled_raw:      {
%            .maxsum_tstart: (optional 1xncells offset in seconds)
%            .force_tmax_optim: dummy field.  if present will trigger use of meea_maxsum_tstart in process (field will be removed afterwards)
% outputs: analysis_out structure variable will be a copy of analysis in with the following fields added
%            .trange:            }
%            .resample_rate:     }} if not present will be set by defaults of 'meea_resample'
%            .resampled_raw:     }
%            .maxsum_tmax.raw (standard) = {nvalves}(n_repeats, ncells) with times chosen by maxsum_monotonic
%                        .filtered (if .resample_filtered is present)
%            .maxsum_rmax.raw (standard) = {nvalves}(n_repeats, ncells) with rates @ tmax chosen by maxsum_monotonic
%                        .filtered (if .resample_filtered is present)
%            *.maxsum_rwindow.raw.s%d (standard) = {nvalves}(n_repeats, ncells) with rates @ tmax chosen by maxsum_monotonic
%                         .filtered.s%d (if .resample_filtered is present)
%                 * "s%d" notation indicates the integration window around
%                         the average tmax for the cell.  Default 's1' for
%                         1 second integration windows,. (CURRENTLY 1 SEC HARD-CODED ONLY)
%            .maxsum_rcumulative = {nvalves}(n_repeats, ncells)
%                         --> contains mean(firingrate(maxsum_tstart:maxsumtmax)) info
%                         
%
% See also multi_electrode_ephys_analyze, ephys, maxsum_monotonic

% Copyright 2009 Julian P. Meeks (Timothy Holy Laboratory)
% 
%% Version control
version = 1.099;
if isfield(analysis_in, 'version')
    if isfield(analysis_in.version, 'maxsum')
        if analysis_in.version.maxsum >= version
            fprintf('meea_maxsum is up to date.\n');
            analysis_out = analysis_in;
            return;
        else
            fprintf('meea_maxsum is out of date. Updating from v. %.2d to v. %.2d.\n', analysis_in.version.maxsum, version);
        end
    end
end

%% Check for variables
analysis_out = analysis_in;

if ~isfield(analysis_out, 'trange') || ~isfield(analysis_out, 'resample_rate') || ~isfield(analysis_out, 'resampled_raw')
    analysis_out = meea_resample(ephysin, analysis_out);
end
if isfield(analysis_out, 'force_deltar_option') % placeholder struct field to indicate choice (for batch processing ease)
    if analysis_out.force_deltar_option
        analysis_out = meea_baserate(ephysin, analysis_out);
    end
    analysis_out = rmfield(analysis_out, 'force_deltar_option');
end
if ~isfield(analysis_out, 'maxsum_valvegroups') || ~isfield(analysis_out, 'maxsum_valveids')
    analysis_out = meea_maxsum_valvegroup(ephysin, analysis_out);
end
if isfield(analysis_out, 'force_tmax_optim') % placeholder struct field to indicate choice (for batch processing ease)
    if analysis_out.force_tmax_optim
        analysis_out = meea_maxsum_tstart(ephysin, analysis_out);
    end
    analysis_out = rmfield(analysis_out, 'force_tmax_optim');
end

%% Fxn main
toffset = analysis_out.trange(1);
resample = analysis_out.resample_rate;

ncells = size(analysis_out.resampled_raw{1},3);
nvalves = size(analysis_out.resampled_raw,2);
ngroups = size(analysis_out.maxsum_valvegroups,2);

if isfield(analysis_out, 'maxsum_tstart')
    startsamp = abs(toffset*resample)+analysis_out.maxsum_tstart*resample; % modify per cell based on responses
    if strmatch('filtered', analysis_out.maxsum_tstart_basis)
        filter_on = true;
    else
        filter_on = false;
    end
else
    startsamp(1:ncells) = abs(toffset*resample);
    filter_on = false;
end

if isfield(analysis_out, 'resampled_deltar')
    deltar_on = true;
else
    deltar_on = false;
end

for ci = 1:ncells
    for gi = 1:ngroups
        totalt = size(analysis_out.resampled_raw{gi},2);
        evalt = size(analysis_out.resampled_raw{gi}(:,startsamp(ci):end,ci),2)-1;
        valvegroups = analysis_out.maxsum_valvegroups{gi};
        ntrials = size(analysis_out.resampled_raw{valvegroups(1)},1);
        for ti = 1:ntrials
            nconcs = size(valvegroups,2);
            this_fr_mean = mean(cell2mat(analysis_out.resampled_raw(valvegroups)),1);
            this_fr = cell2mat(analysis_out.resampled_raw(valvegroups));
            if filter_on == true
                this_filt_fr_mean = mean(cell2mat(analysis_out.resampled_filter(valvegroups)),1);
                this_filt_fr = cell2mat(analysis_out.resampled_filter(valvegroups));
            end
            if deltar_on == true
                this_dr_mean = mean(cell2mat(analysis_out.resampled_deltar(valvegroups)),1);
                this_dr = cell2mat(analysis_out.resampled_deltar(valvegroups));
                if filter_on == true
                    this_filt_dr_mean = mean(cell2mat(analysis_out.resampled_deltar_filtered(valvegroups)),1);
                    this_filt_dr = cell2mat(analysis_out.resampled_deltar_filtered(valvegroups));
                end
            end
            if size(this_fr,2) > totalt
                temp_fr = this_fr;
                temp_fr_mean = this_fr_mean;
                this_fr = [];
                this_fr_mean = [];
                if filter_on == true
                    temp_filt_fr = this_filt_fr;
                    temp_filt_fr_mean = this_filt_fr_mean;
                    this_filt_fr = [];
                    this_filt_fr_mean = [];
                end
                if deltar_on == true
                    temp_dr = this_dr;
                    temp_dr_mean = this_dr_mean;
                    this_dr = [];
                    this_dr_mean = [];
                    if filter_on == true
                        temp_filt_dr = this_filt_dr;
                        temp_filt_dr_mean = this_filt_dr_mean;
                        this_filt_dr = [];
                        this_filt_dr_mean = [];
                    end
                end
                count = 1;
                while ~isempty(temp_fr)
                    this_fr(count,:) = temp_fr(ti,1:totalt,ci);
                    this_fr_mean(count,:) = temp_fr_mean(:,1:totalt,ci);
                    temp_fr(:,1:totalt,:) = [];
                    temp_fr_mean(:,1:totalt,:) = [];
                    if filter_on == true
                        this_filt_fr(count,:) = temp_filt_fr(ti,1:totalt,ci);
                        this_filt_fr_mean(count,:) = temp_filt_fr_mean(:,1:totalt,ci);
                        temp_filt_fr(:,1:totalt,:) = [];
                        temp_filt_fr_mean(:,1:totalt,:) = [];
                    end
                    if deltar_on == true
                        this_dr(count,:) = temp_dr(ti,1:totalt,ci);
                        this_dr_mean(count,:) = temp_dr_mean(:,1:totalt,ci);
                        temp_dr(:,1:totalt,:) = [];
                        temp_dr_mean(:,1:totalt,:) = [];
                        if filter_on == true
                            this_filt_dr(count,:) = temp_filt_dr(ti,1:totalt,ci);
                            this_filt_dr_mean(count,:) = temp_filt_dr_mean(:,1:totalt,ci);
                            temp_filt_dr(:,1:totalt,:) = [];
                            temp_filt_dr_mean(:,1:totalt,:) = [];
                        end
                    end
                    count = count+1;
                end
            else
                this_fr = this_fr(ti,:,ci);
                this_fr_mean = this_fr_mean(:,:,ci);
                if filter_on == true
                    this_filt_fr = this_filt_fr(ti,:,ci);
                    this_filt_fr_mean = this_filt_fr_mean(:,:,ci);
                end
            end            
            fr_t = maxsum_monotonic(this_fr_mean(:,startsamp(ci):end)');
            if filter_on == true
                fr_filt_t = maxsum_monotonic(this_filt_fr_mean(:,startsamp(ci):end)');
            end
            %% note: if you wonder where maxsum_monotonic for
            %% deltar-based set is: the fr_tmax and fr_filt_tmax have
            %% same values as dr_tmax, dr_filt_tmax
            % finally assign values
            for conc_idx = 1:nconcs
                fr_tmax{valvegroups(conc_idx)}(ti,ci) = fr_t(conc_idx);
                fr_rmax{valvegroups(conc_idx)}(ti,ci) = analysis_out.resampled_raw{valvegroups(conc_idx)}(ti,startsamp(ci)+fr_t(conc_idx)-1,ci);                
                fr_samplerange.s1 = startsamp(ci)+fr_t(conc_idx)-(resample):startsamp(ci)+fr_t(conc_idx)-1;
                if fr_t(conc_idx) <= resample
                    fr_samplerange.cumul = startsamp(ci)+fr_t(conc_idx)-(resample):startsamp(ci)+fr_t(conc_idx)-1;
                else
                    fr_samplerange.cumul = startsamp(ci):startsamp(ci)+fr_t(conc_idx)-1;
                end
                fr_rwindow.s1{valvegroups(conc_idx)}(ti,ci) = mean(analysis_out.resampled_raw{valvegroups(conc_idx)}(ti,fr_samplerange.s1, ci));% hard-coded for now, can be made an option.  .s1 = 1sec integration;
                fr_rcumulative{valvegroups(conc_idx)}(ti,ci) = mean(analysis_out.resampled_raw{valvegroups(conc_idx)}(ti,fr_samplerange.cumul, ci));
                if filter_on == true
                    fr_filt_tmax{valvegroups(conc_idx)}(ti,ci) = fr_filt_t(conc_idx);
                    fr_filt_rmax{valvegroups(conc_idx)}(ti,ci) = analysis_out.resampled_filter{valvegroups(conc_idx)}(ti,startsamp(ci)+fr_filt_t(conc_idx)-1,ci);
                    filt_fr_samplerange.s1 = startsamp(ci)+fr_filt_t(conc_idx)-(resample):startsamp(ci)+fr_filt_t(conc_idx)-1;  % startsamp(ci)+fr_filt_t(conc_idx)-floor(resample/2):startsamp(ci)+fr_filt_t(conc_idx)+floor(resample/2-1);
                    if fr_filt_t(conc_idx) <= resample
                        filt_fr_samplerange.cumul = startsamp(ci)+fr_filt_t(conc_idx)-(resample):startsamp(ci)+fr_filt_t(conc_idx)-1;
                    else
                        filt_fr_samplerange.cumul = startsamp(ci):startsamp(ci)+fr_filt_t(conc_idx)-1;
                    end
                    fr_filt_rwindow.s1{valvegroups(conc_idx)}(ti,ci) = mean(analysis_out.resampled_filter{valvegroups(conc_idx)}(ti,filt_fr_samplerange.s1, ci));
                    fr_filt_rcumulative{valvegroups(conc_idx)}(ti,ci) = mean(analysis_out.resampled_filter{valvegroups(conc_idx)}(ti,fr_samplerange.cumul, ci));
                end
                if deltar_on == true
                    dr_tmax{valvegroups(conc_idx)}(ti,ci) = fr_t(conc_idx);
                    dr_rmax{valvegroups(conc_idx)}(ti,ci) = analysis_out.resampled_deltar{valvegroups(conc_idx)}(ti,startsamp(ci)+fr_t(conc_idx)-1,ci);
                    dr_rwindow.s1{valvegroups(conc_idx)}(ti,ci) = mean(analysis_out.resampled_deltar{valvegroups(conc_idx)}(ti,fr_samplerange.s1, ci));
                    dr_rcumulative{valvegroups(conc_idx)}(ti,ci) = mean(analysis_out.resampled_deltar{valvegroups(conc_idx)}(ti,fr_samplerange.cumul, ci));
                    if filter_on == true
                        dr_filt_tmax{valvegroups(conc_idx)}(ti,ci) = fr_filt_t(conc_idx);
                        dr_filt_rmax{valvegroups(conc_idx)}(ti,ci) = analysis_out.resampled_deltar_filtered{valvegroups(conc_idx)}(ti,startsamp(ci)+fr_filt_t(conc_idx)-1,ci);
                        dr_filt_rwindow.s1{valvegroups(conc_idx)}(ti,ci) = mean(analysis_out.resampled_deltar_filtered{valvegroups(conc_idx)}(ti,filt_fr_samplerange.s1, ci));
                        dr_filt_rcumulative{valvegroups(conc_idx)}(ti,ci) = mean(analysis_out.resampled_deltar_filtered{valvegroups(conc_idx)}(ti,filt_fr_samplerange.cumul, ci));
                    end
                end
            end
        end
    end
end
%% Compile means, sems, ranksum, ttest values:

ringvalve = strmatch('Ringer''s',analysis_out.valvelabels);
if isempty(ringvalve)
    error('No Ringer''s valve found');
end

for vi = 1:nvalves
    for ci = 1:ncells
        %% means
        maxsum.mean.fr.raw(vi,ci) = mean(fr_rmax{vi}(:,ci),1);
        maxsum.mean.frwindow.raw.s1(vi,ci) = mean(fr_rwindow.s1{vi}(:,ci),1);
        maxsum.mean.frcumul.raw(vi,ci) = mean(fr_rcumulative{vi}(:,ci),1);
        if filter_on == true
            maxsum.mean.fr.filtered(vi,ci) = mean(fr_filt_rmax{vi}(:,ci),1);
            maxsum.mean.frwindow.filtered.s1(vi,ci) = mean(fr_filt_rwindow.s1{vi}(:,ci),1);
            maxsum.mean.frcumul.filtered(vi,ci) = mean(fr_filt_rcumulative{vi}(:,ci),1);
        end
        if deltar_on == true
            maxsum.mean.dr.raw(vi,ci) = mean(dr_rmax{vi}(:,ci),1);
            maxsum.mean.drwindow.raw.s1(vi,ci) = mean(dr_rwindow.s1{vi}(:,ci),1);
            maxsum.mean.drcumul.raw(vi,ci) = mean(dr_rcumulative{vi}(:,ci),1);
            if filter_on == true
                maxsum.mean.dr.filtered(vi,ci) = mean(dr_filt_rmax{vi}(:,ci),1);
                maxsum.mean.drwindow.filtered.s1(vi,ci) = mean(dr_filt_rwindow.s1{vi}(:,ci),1);
                maxsum.mean.drcumul.filtered(vi,ci) = mean(dr_filt_rcumulative{vi}(:,ci),1);
            end
        end
        %% sems
        maxsum.sem.fr.raw(vi,ci) = std(fr_rmax{vi}(:,ci),0,1)/sqrt(size(fr_rmax{vi}(:,ci),1));
        maxsum.sem.frwindow.raw.s1(vi,ci) = std(fr_rwindow.s1{vi}(:,ci),0,1)/sqrt(size(fr_rwindow.s1{vi}(:,ci),1));
        maxsum.sem.frcumul.raw(vi,ci) = std(fr_rcumulative{vi}(:,ci),0,1)/sqrt(size(fr_rcumulative{vi}(:,ci),1));
        if filter_on == true
            maxsum.sem.fr.filtered(vi,ci) = std(fr_filt_rmax{vi}(:,ci),0,1)/sqrt(size(fr_filt_rmax{vi}(:,ci),1));
            maxsum.sem.frwindow.filtered.s1(vi,ci) = std(fr_filt_rwindow.s1{vi}(:,ci),0,1)/sqrt(size(fr_filt_rwindow.s1{vi}(:,ci),1));
            maxsum.sem.frcumul.filtered(vi,ci) = std(fr_filt_rcumulative{vi}(:,ci),0,1)/sqrt(size(fr_filt_rcumulative{vi}(:,ci),1));
        end
        if deltar_on == true
            maxsum.sem.dr.raw(vi,ci) = std(dr_rmax{vi}(:,ci),0,1)/sqrt(size(dr_rmax{vi}(:,ci),1));
            maxsum.sem.drwindow.raw.s1(vi,ci) = std(dr_rwindow.s1{vi}(:,ci),0,1)/sqrt(size(dr_rwindow.s1{vi}(:,ci),1));
            maxsum.sem.drcumul.raw(vi,ci) = std(dr_rcumulative{vi}(:,ci),0,1)/sqrt(size(dr_rcumulative{vi}(:,ci),1));
            if filter_on == true
                maxsum.sem.dr.filtered(vi,ci) = std(dr_filt_rmax{vi}(:,ci),0,1)/sqrt(size(dr_filt_rmax{vi}(:,ci),1));
                maxsum.sem.drwindow.filtered.s1(vi,ci) = std(dr_filt_rwindow.s1{vi}(:,ci),0,1)/sqrt(size(dr_filt_rwindow.s1{vi}(:,ci),1));
                maxsum.sem.drcumul.filtered(vi,ci) = std(dr_filt_rcumulative{vi}(:,ci),0,1)/sqrt(size(dr_filt_rcumulative{vi}(:,ci),1));
            end
        end
        %% ranksum
        maxsum.ranksum.fr.raw(vi,ci) = ranksum(fr_rmax{vi}(:,ci),fr_rmax{ringvalve}(:,ci));
        maxsum.ranksum.frwindow.raw.s1(vi,ci) = ranksum(fr_rwindow.s1{vi}(:,ci),fr_rwindow.s1{ringvalve}(:,ci));
        maxsum.ranksum.frcumul.raw(vi,ci) = ranksum(fr_rcumulative{vi}(:,ci),fr_rcumulative{ringvalve}(:,ci));
        if filter_on == true
            maxsum.ranksum.fr.filtered(vi,ci) = ranksum(fr_filt_rmax{vi}(:,ci),fr_filt_rmax{ringvalve}(:,ci));
            maxsum.ranksum.frwindow.filtered.s1(vi,ci) = ranksum(fr_filt_rwindow.s1{vi}(:,ci),fr_filt_rwindow.s1{ringvalve}(:,ci));
            maxsum.ranksum.frcumul.filtered(vi,ci) = ranksum(fr_filt_rcumulative{vi}(:,ci),fr_filt_rcumulative{ringvalve}(:,ci));
        end
        if deltar_on == true
            maxsum.ranksum.dr.raw(vi,ci) = ranksum(dr_rmax{vi}(:,ci),dr_rmax{ringvalve}(:,ci));
            maxsum.ranksum.drwindow.raw.s1(vi,ci) = ranksum(dr_rwindow.s1{vi}(:,ci),dr_rwindow.s1{ringvalve}(:,ci));
            maxsum.ranksum.drcumul.raw(vi,ci) = ranksum(dr_rcumulative{vi}(:,ci),dr_rcumulative{ringvalve}(:,ci));
            if filter_on == true
                maxsum.ranksum.dr.filtered(vi,ci) = ranksum(dr_filt_rmax{vi}(:,ci),dr_filt_rmax{ringvalve}(:,ci));
                maxsum.ranksum.drwindow.filtered.s1(vi,ci) = ranksum(dr_filt_rwindow.s1{vi}(:,ci),dr_filt_rwindow.s1{ringvalve}(:,ci));
                maxsum.ranksum.drcumul.filtered(vi,ci) = ranksum(dr_filt_rcumulative{vi}(:,ci),dr_filt_rcumulative{ringvalve}(:,ci));
            end
        end
        %% ttest p
        [temp maxsum.ttest.fr.raw(vi,ci)] = ttest2(fr_rmax{vi}(:,ci),fr_rmax{ringvalve}(:,ci));
        [temp maxsum.ttest.frwindow.raw.s1(vi,ci)] = ttest2(fr_rwindow.s1{vi}(:,ci),fr_rwindow.s1{ringvalve}(:,ci));
        [temp maxsum.ttest.frcumul.raw(vi,ci)] = ttest2(fr_rcumulative{vi}(:,ci),fr_rcumulative{ringvalve}(:,ci));
        if filter_on == true
            [temp maxsum.ttest.fr.filtered(vi,ci)] = ttest2(fr_filt_rmax{vi}(:,ci),fr_filt_rmax{ringvalve}(:,ci));
            [temp maxsum.ttest.frwindow.filtered.s1(vi,ci)] = ttest2(fr_filt_rwindow.s1{vi}(:,ci),fr_filt_rwindow.s1{ringvalve}(:,ci));
            [temp maxsum.ttest.frcumul.filtered(vi,ci)] = ttest2(fr_filt_rcumulative{vi}(:,ci),fr_filt_rcumulative{ringvalve}(:,ci));
        end
        if deltar_on == true
            [temp maxsum.ttest.dr.raw(vi,ci)] = ttest2(dr_rmax{vi}(:,ci),dr_rmax{ringvalve}(:,ci));
            [temp maxsum.ttest.drwindow.raw.s1(vi,ci)] = ttest2(dr_rwindow.s1{vi}(:,ci),dr_rwindow.s1{ringvalve}(:,ci));
            [temp maxsum.ttest.drcumul.raw(vi,ci)] = ttest2(dr_rcumulative{vi}(:,ci),dr_rcumulative{ringvalve}(:,ci));
            if filter_on == true
                [temp maxsum.ttest.dr.filtered(vi,ci)] = ttest2(dr_filt_rmax{vi}(:,ci),dr_filt_rmax{ringvalve}(:,ci));
                [temp maxsum.ttest.drwindow.filtered.s1(vi,ci)] = ttest2(dr_filt_rwindow.s1{vi}(:,ci),dr_filt_rwindow.s1{ringvalve}(:,ci));
                [temp maxsum.ttest.drcumul.filtered(vi,ci)] = ttest2(dr_filt_rcumulative{vi}(:,ci),dr_filt_rcumulative{ringvalve}(:,ci));
            end
        end
        %% means
        maxsum.mean.tmax.raw(vi,ci) = mean(fr_tmax{vi}(:,ci),1);
        if filter_on == true
            maxsum.mean.tmax.filtered(vi,ci) = mean(fr_filt_tmax{vi}(:,ci),1);
        end
    end
end

%% Finish assignments
maxsum.values.fr.raw = fr_rmax;
maxsum.values.frwindow.raw = fr_rwindow;
maxsum.values.frcumul.raw = fr_rcumulative;
maxsum.values.tmax.raw = fr_tmax;
if filter_on == true
    maxsum.values.fr.filtered = fr_filt_rmax;
    maxsum.values.frwindow.filtered = fr_filt_rwindow;
    maxsum.values.frcumul.filtered = fr_filt_rcumulative;
    maxsum.values.tmax.filtered = fr_filt_tmax;
end
if deltar_on == true
    maxsum.values.dr.raw = dr_rmax;
    maxsum.values.drwindow.raw = dr_rwindow;
    maxsum.values.drcumul.raw = dr_rcumulative;
    if filter_on == true
        maxsum.values.dr.filtered = dr_filt_rmax;
        maxsum.values.drwindow.filtered = dr_filt_rwindow;
        maxsum.values.drcumul.filtered = dr_filt_rcumulative;
    end
end


analysis_out.maxsum = maxsum;

analysis_out.version.maxsum = version;

end