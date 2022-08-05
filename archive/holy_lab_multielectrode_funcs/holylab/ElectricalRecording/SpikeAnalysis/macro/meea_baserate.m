function analysis_out = meea_baserate(ephysin, analysis_in)
% meea_baserate calculates baseline firing rates for each cell in input
% Syntax: analysis_out = meea_baserate(ephysin, analysis_in)
%
% Inputs:  ephysin: multi electrode ephys structure
%          analysis_in: analysis structure with required subfields:
%            trange:             {
%            resample_rate:     {{ if not present, will be set to defaults of 'meea_resample'
%            resampled_raw:      {
%            
% Outputs: analysis_out structure variable will be a copy of analysis in with the following fields added
%            .trange:            }
%            .resample_rate:     }} if not present will be set by defaults of 'meea_resample'
%            .resampled_raw      }
%            .resample_filter:   } (only if .resampled_filter subfield present in analysis_in)
%            .resampled_deltar {nvalves}(nrepeats, nsamples, ncells) with baseline subtracted
%            .resampled_deltar_filtered (if resample_filter present)
%            .baserate: (n_cells) scalar firing rate during trange(1):0 across ALL TRIALS
%            .baserate_filtered: (if resample_filter present)
%
% See also multi_electrode_ephys_analyze, meea_resample

% Copyright 2009 Julian P. Meeks (Timothy Holy Laboratory)
% 
% Version History:
% 2009_01-2009_05: Wrote it (JPM)
% 2009_08_13: checked and fixed help info (JPM)
%% Version control
version = 1.01;
if isfield(analysis_in, 'version')
    if isfield(analysis_in.version, 'baserate')
        if analysis_in.version.baserate == version
            fprintf('meea_baserate is up to date.\n');
            analysis_out = analysis_in;
            return;
        else
            fprintf('meea_baserate is out of date. Updating from v. %2d to v. %2d.\n', analysis_in.version.baserate, version);
        end
    end
end

%% Fxn main:
analysis_out = analysis_in;

if ~isfield(analysis_out, 'trange') || ~isfield(analysis_out, 'resample_rate') || ~isfield(analysis_out, 'resampled_raw')
    analysis_out = meea_resample(ephysin, analysis_out);
end

if isfield(analysis_out, 'resampled_filter')
    filter_on = true;
else
    filter_on = false;
end

resample = analysis_out.resample_rate;
toffset = analysis_out.trange(1);
ncells = size(analysis_out.resampled_raw{1},3);
nvalves = size(analysis_out.resampled_raw,2);
for ci = 1:ncells
    for vi = 1:nvalves
        baseline_matrix(ci,vi,:) = mean(analysis_out.resampled_raw{vi}(:,1:abs(toffset)*resample,ci),1);
        if filter_on == true
            filtered_baseline_matrix(ci,vi,:) = mean(analysis_out.resampled_filter{vi}(:,1:abs(toffset)*resample, ci),1);
        end
    end
end
avg_baseline_matrix = mean(mean(baseline_matrix,3),2);
if filter_on == true
    avg_filtered_baseline_matrix = mean(mean(filtered_baseline_matrix,3),2);
end
for vi = 1:nvalves
    for ci = 1:ncells
        resampled_deltar{vi}(:,:,ci) = analysis_out.resampled_raw{vi}(:,:,ci)-avg_baseline_matrix(ci);
        if filter_on == true
            filtered_resampled_deltar{vi}(:,:,ci) = analysis_out.resampled_filter{vi}(:,:,ci)-avg_filtered_baseline_matrix(ci);
        end
    end
end


%% Finish assignments
analysis_out.baserate = avg_baseline_matrix';
analysis_out.resampled_deltar = resampled_deltar;
if filter_on == true;
   analysis_out.baserate_filtered = avg_filtered_baseline_matrix';
   analysis_out.resampled_deltar_filtered = filtered_resampled_deltar;
end

analysis_out.version.baserate = version;

end