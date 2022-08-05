function analysis_out = meea_maxsum_timeweight(ephysin, analysis_in)
% meea_maxsum_timeweight modifies deltar_max vals based on early/late tmax
%
% **NOTE: this analysis has fallen out of favor, consider meea_maxsum:
%        --> analysis.maxsum.mean.drcumulative.(raw or filtered)
%
% syntax: analysis_out = meea_maxsum_timeweight(ephysin, analysis_in)
%
% inputs:  ephysin: multi electrode ephys structure
%          analysis_in: analysis structure with optional subfields:
%            trange:             {
%            resample_rate:     {{ if not present, will be set to defaults of 'meea_resample'
%            resampled_raw:      {
%            maxsum...          {{ if not present, will be calculated by meea_maxsum
% outputs: analysis_out structure variable will be a copy of analysis in with the following fields added
%            (*all values assigned by precursor functions, see respective functions for subfield definitions
%            .maxsum_timeweight.values.. all subfields are {valves}(trials,cells)
%                              .fr.raw (uses raw firing rate and tmax on raw firing rate)
%                              .fr.filtered (uses filtered firing rate and tmax on filtered firing rate)
%                              .dr.raw (uses deltar and tmax on raw firing rate)
%                              .dr.filtered (uses deltar filtered and tmax on filtered firing rate)
%            .maxsum_timeweight.mean (mean over trials using same subfields as above) (nvalves, ncells)
%            .maxsum_timeweight.sem (standard error using same subfields above) (nvalves, ncells)
%            .maxsum_timeweight.ranksum (ranksum rho using same subfields above) (nvalves, ncells)
%            .maxsum_timeweight.ttest (ttest p value using same subfields above) (nvalves, ncells)
%
% See also multi_electrode_ephys_analyze, ephys, meea_resample, meea_maxsum, meea_baseline, meea_maxsum_tstart

% Copyright 2009 Julian P. Meeks (Timothy Holy Laboratory)
%
% Version History:
% 2009_01-2009_05: Wrote it (JPM)
% 2009_08_14: checked and fixed help info (JPM) 
%% Version control
version = 1.099;
if isfield(analysis_in, 'version')
    if isfield(analysis_in.version, 'maxsum_timeweight')
        if analysis_in.version.maxsum_timeweight >= version
            fprintf('meea_maxsum_timeweight is up to date.\n');
            analysis_out = analysis_in;
            return;
        else
            fprintf('meea_maxsum_timeweight is out of date. Updating from v. %.2d to v. %.2d.\n', analysis_in.version.maxsum_tstart, version);
        end
    end
end

%% Check for variables
analysis_out = analysis_in;

if ~isfield(analysis_out, 'maxsum') || analysis_in.version.maxsum < version
    analysis_out = meea_maxsum(ephysin, analysis_out);
end

%% Fxn main
toffset = analysis_out.trange(1);
resample = analysis_out.resample_rate;
ncells = size(analysis_out.resampled_raw{1},3);
nvalves = size(analysis_out.resampled_raw,2);
totalt = size(analysis_out.resampled_raw{1},2)-abs(toffset)*resample;
startsamp = analysis_out.maxsum_tstart*resample;
evalt = totalt-startsamp;

if isfield(analysis_out.maxsum.values, 'dr')
    deltar_on = true;
else
    deltar_on = false;
end
if isfield(analysis_out.maxsum.mean.fr, 'filtered')
    filter_on = true;
else
    filter_on = false;
end

for vi = 1:nvalves
    ntrials = size(analysis_out.resampled_raw{vi},1);
    for ci = 1:ncells
        for ti = 1:ntrials
            timeweight.values.fr.raw{vi}(ti,ci) = analysis_out.maxsum.values.fr.raw{vi}(ti,ci)*(1-analysis_out.maxsum.values.tmax.raw{vi}(ti,ci)/evalt(ci));
            timeweight.values.frwindow.raw.s1{vi}(ti,ci) = analysis_out.maxsum.values.frwindow.raw.s1{vi}(ti,ci)*(1-analysis_out.maxsum.values.tmax.raw{vi}(ti,ci)/evalt(ci));
            if filter_on == true
                timeweight.values.fr.filtered{vi}(ti,ci) = analysis_out.maxsum.values.fr.filtered{vi}(ti,ci)*(1-analysis_out.maxsum.values.tmax.filtered{vi}(ti,ci)/evalt(ci));
                timeweight.values.frwindow.filtered.s1{vi}(ti,ci) = analysis_out.maxsum.values.frwindow.filtered.s1{vi}(ti,ci)*(1-analysis_out.maxsum.values.tmax.filtered{vi}(ti,ci)/evalt(ci));
            end
            if deltar_on == true
                timeweight.values.dr.raw{vi}(ti,ci) = analysis_out.maxsum.values.dr.raw{vi}(ti,ci)*(1-analysis_out.maxsum.values.tmax.raw{vi}(ti,ci)/evalt(ci));
                timeweight.values.drwindow.raw.s1{vi}(ti,ci) = analysis_out.maxsum.values.drwindow.raw.s1{vi}(ti,ci)*(1-analysis_out.maxsum.values.tmax.raw{vi}(ti,ci)/evalt(ci));
                if filter_on == true
                    timeweight.values.dr.filtered{vi}(ti,ci) = analysis_out.maxsum.values.dr.filtered{vi}(ti,ci)*(1-analysis_out.maxsum.values.tmax.filtered{vi}(ti,ci)/evalt(ci));
                    timeweight.values.drwindow.filtered.s1{vi}(ti,ci) = analysis_out.maxsum.values.drwindow.filtered.s1{vi}(ti,ci)*(1-analysis_out.maxsum.values.tmax.filtered{vi}(ti,ci)/evalt(ci));
                end
            end
        end
        %% means
        timeweight.mean.fr.raw(vi,ci) = mean(timeweight.values.fr.raw{vi}(:,ci),1);
        timeweight.mean.frwindow.raw.s1(vi,ci) = mean(timeweight.values.frwindow.raw.s1{vi}(:,ci),1);
        if filter_on == true
            timeweight.mean.fr.filtered(vi,ci) = mean(timeweight.values.fr.filtered{vi}(:,ci),1);
            timeweight.mean.frwindow.filtered.s1(vi,ci) = mean(timeweight.values.frwindow.filtered.s1{vi}(:,ci),1);
        end
        if deltar_on == true
            timeweight.mean.dr.raw(vi,ci) = mean(timeweight.values.dr.raw{vi}(:,ci),1);
            timeweight.mean.drwindow.raw.s1(vi,ci) = mean(timeweight.values.drwindow.raw.s1{vi}(:,ci),1);
            if filter_on == true
                timeweight.mean.dr.filtered(vi,ci) = mean(timeweight.values.dr.filtered{vi}(:,ci),1);
                timeweight.mean.drwindow.filtered.s1(vi,ci) = mean(timeweight.values.drwindow.filtered.s1{vi}(:,ci),1);
            end
        end
        %% sems
        timeweight.sem.fr.raw(vi,ci) = std(timeweight.values.fr.raw{vi}(:,ci),0,1)/sqrt(ntrials);
        timeweight.sem.frwindow.raw.s1(vi,ci) = std(timeweight.values.frwindow.raw.s1{vi}(:,ci),0,1)/sqrt(ntrials);
        if filter_on == true
            timeweight.sem.fr.filtered(vi,ci) = std(timeweight.values.fr.filtered{vi}(:,ci),0,1)/sqrt(ntrials);
            timeweight.sem.frwindow.filtered.s1(vi,ci) = std(timeweight.values.frwindow.filtered.s1{vi}(:,ci),0,1)/sqrt(ntrials);
        end
        if deltar_on == true
            timeweight.sem.dr.raw(vi,ci) = std(timeweight.values.dr.raw{vi}(:,ci),0,1)/sqrt(ntrials);
            timeweight.sem.drwindow.raw.s1(vi,ci) = std(timeweight.values.drwindow.raw.s1{vi}(:,ci),0,1)/sqrt(ntrials);
            if filter_on == true
                timeweight.sem.dr.filtered(vi,ci) = std(timeweight.values.dr.filtered{vi}(:,ci),0,1)/sqrt(ntrials);
                timeweight.sem.drwindow.filtered.s1(vi,ci) = std(timeweight.values.drwindow.filtered.s1{vi}(:,ci),0,1)/sqrt(ntrials);
            end
        end
    end
end

ringvalve = strmatch('Ringer''s',analysis_out.valvelabels);
if isempty(ringvalve)
    error('No Ringer''s valve found');
end

% do stats
for vi = 1:nvalves
    for ci = 1:ncells
        %% ranksum
        timeweight.ranksum.fr.raw(vi,ci) = ranksum(timeweight.values.fr.raw{vi}(:,ci),timeweight.values.fr.raw{ringvalve}(:,ci));
        timeweight.ranksum.frwindow.raw.s1(vi,ci) = ranksum(timeweight.values.frwindow.raw.s1{vi}(:,ci),timeweight.values.frwindow.raw.s1{ringvalve}(:,ci));
        if filter_on == true
            timeweight.ranksum.fr.filtered(vi,ci) = ranksum(timeweight.values.fr.filtered{vi}(:,ci),timeweight.values.fr.filtered{ringvalve}(:,ci));
            timeweight.ranksum.frwindow.filtered.s1(vi,ci) = ranksum(timeweight.values.frwindow.filtered.s1{vi}(:,ci),timeweight.values.frwindow.filtered.s1{ringvalve}(:,ci));
        end
        if deltar_on == true
            timeweight.ranksum.dr.raw(vi,ci) = ranksum(timeweight.values.dr.raw{vi}(:,ci),timeweight.values.dr.raw{ringvalve}(:,ci));
            timeweight.ranksum.drwindow.raw.s1(vi,ci) = ranksum(timeweight.values.drwindow.raw.s1{vi}(:,ci),timeweight.values.drwindow.raw.s1{ringvalve}(:,ci));
            if filter_on == true
                timeweight.ranksum.dr.filtered(vi,ci) = ranksum(timeweight.values.dr.filtered{vi}(:,ci),timeweight.values.dr.filtered{ringvalve}(:,ci));
                timeweight.ranksum.drwindow.filtered.s1(vi,ci) = ranksum(timeweight.values.drwindow.filtered.s1{vi}(:,ci),timeweight.values.drwindow.filtered.s1{ringvalve}(:,ci));
            end
        end
        %% ttest p
        [temp timeweight.ttest.fr.raw(vi,ci)] = ttest2(timeweight.values.fr.raw{vi}(:,ci),timeweight.values.fr.raw{ringvalve}(:,ci));
        [temp timeweight.ttest.frwindow.raw.s1(vi,ci)] = ttest2(timeweight.values.frwindow.raw.s1{vi}(:,ci),timeweight.values.frwindow.raw.s1{ringvalve}(:,ci));
        if filter_on == true
            [temp timeweight.ttest.fr.filtered(vi,ci)] = ttest2(timeweight.values.fr.filtered{vi}(:,ci),timeweight.values.fr.filtered{ringvalve}(:,ci));
            [temp timeweight.ttest.frwindow.filtered.s1(vi,ci)] = ttest2(timeweight.values.frwindow.filtered.s1{vi}(:,ci),timeweight.values.frwindow.filtered.s1{ringvalve}(:,ci));
        end
        if deltar_on == true
            [temp timeweight.ttest.dr.raw(vi,ci)] = ttest2(timeweight.values.dr.raw{vi}(:,ci),timeweight.values.dr.raw{ringvalve}(:,ci));
            [temp timeweight.ttest.drwindow.raw.s1(vi,ci)] = ttest2(timeweight.values.drwindow.raw.s1{vi}(:,ci),timeweight.values.drwindow.raw.s1{ringvalve}(:,ci));
            if filter_on == true
                [temp timeweight.ttest.dr.filtered(vi,ci)] = ttest2(timeweight.values.dr.filtered{vi}(:,ci),timeweight.values.dr.filtered{ringvalve}(:,ci));
                [temp timeweight.ttest.drwindow.filtered.s1(vi,ci)] = ttest2(timeweight.values.drwindow.filtered.s1{vi}(:,ci),timeweight.values.drwindow.filtered.s1{ringvalve}(:,ci));
            end
        end
    end
end

%% Finish assignments
analysis_out.maxsum_timeweight = timeweight;

analysis_out.version.maxsum_timeweight = version;

end