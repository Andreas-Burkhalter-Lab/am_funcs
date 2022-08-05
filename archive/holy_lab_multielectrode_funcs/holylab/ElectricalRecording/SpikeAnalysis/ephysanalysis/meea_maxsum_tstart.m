function analysis_out = meea_maxsum_tstart(ephysin, analysis_in)
% meea_maxsum_tstart sets time offsets for multielectrode ephys analysis
%
% syntax: analysis_out = meea_maxsum_tstart(ephysin, analysis_in)
%
% inputs:  ephysin: multi electrode ephys structure
%          analysis_in: analysis structure with optional subfields:
%            trange:             {
%            resample_rate:     {{ if not present, will be set to defaults of 'meea_resample'
%            resample_filter:    {
%            resampled_raw:      {
%            resampled_filter:   {
%            min_tstart_offset: gives the minimum integration time prior to tmax (in seconds, default 1s)
% outputs: analysis_out structure variable will be a copy of analysis in with the following fields added
%            .trange:            }
%            .resample_rate:     }} if not present will be set by defaults of 'meea_resample'
%            .resample_filter:   }
%            .resampled_raw:     }
%            .resampled_filter:  } 
%            .min_tstart_offset: default 1s as described above
%            .maxsum_tstart: set by dimensions of
%                            .resampled_filter{{nvalves}(nrepeats,nsamples,ncells)} -->
%                            .maxsum_tstart(ncells x nresample vals)
% See also multi_electrode_ephys_analyze, ephys, cheby1

% Copyright 2009 Julian P. Meeks (Timothy Holy Laboratory)
% 
%% Version control
version = 1.099;
if isfield(analysis_in, 'version')
    if isfield(analysis_in.version, 'maxsum_tstart')
        if analysis_in.version.maxsum_tstart >= version
            fprintf('meea_maxsum_tstart is up to date.\n');
            analysis_out = analysis_in;
            return;
        else
            fprintf('meea_maxsum_tstart is out of date. Updating from v. %2d to v. %2d.\n', analysis_in.version.maxsum_tstart, version);
        end
    end
end

%% just trying this out
tstart = findst_time(ephysin);
filter_on = true;
analysis_out = analysis_in;
% %% Check for variables
% analysis_out = analysis_in;
% 
% if ~isfield(analysis_out, 'trange') || ~isfield(analysis_out, 'resample_rate') || ~isfield(analysis_out, 'resampled_raw')
%     analysis_out = meea_resample(ephysin, analysis_out);
% end
% if ~isfield(analysis_out, 'maxsum_valvegroups') || ~isfield(analysis_out, 'maxsum_valveids')
%     analysis_out = meea_maxsum_valvegroup(ephysin, analysis_out);
% end
% 
% analysis_out = default(analysis_out, 'min_tstart_offset', 1.0);
% 
% %% Fxn main
% toffset = analysis_out.trange(1);
% resample = analysis_out.resample_rate;
% ncells = size(analysis_out.resampled_raw{1},3);
% nvalves = size(analysis_out.resampled_raw,2);
% startsamp = abs(toffset*resample);
% tmax = zeros([nvalves ncells]);
% 
% if isfield(analysis_out, 'resampled_filter')
%     filter_on = true;
% else
%     filter_on = false;
% end
% 
% % set up avg matrix
% for ci = 1:ncells
%     for vi = 1:nvalves
%         if filter_on ==true
%             thisavg(vi,:,ci) = mean(analysis_out.resampled_filter{vi}(:,:,ci),1);
%         else
%             thisavg(vi,:,ci) = mean(analysis_out.resampled_raw{vi}(:,:,ci),1);
%         end
%     end
% end
% 
% n_valvegroups = size(analysis_out.maxsum_valvegroups,2);
% for ci = 1:ncells
%     evalt = size(thisavg(:,startsamp:end,ci),2)-1;
%     for gi = 1:n_valvegroups
%         these_valves = analysis_out.maxsum_valvegroups{gi};
%         tmax(these_valves,ci) = maxsum_monotonic(thisavg(these_valves,startsamp:end,ci)');
%         drmax(these_valves,ci) = diag(thisavg(these_valves,startsamp+tmax(these_valves,ci)-1,ci));
%     end
% end
% 
% for ci = 1:ncells
%     evalt = size(thisavg(:,startsamp:end,ci),2)-1;
%     minresponse(ci) = max([3*std(drmax(:,ci)) 5]);  % 5 Hz deltar is minimum criterion
% %     if max(drmax(:,ci)) < minresponse(ci)
% %         tmax(:,ci) = 0;
% %     else
%         tmax(drmax(:,ci) < minresponse(ci),ci) = evalt;
% %     end
% end
% 
% tstart_per_cell = min(tmax,[],1);
% tstart_per_cell(tstart_per_cell == evalt) = 0;
% tstart_per_cell = tstart_per_cell-analysis_out.min_tstart_offset*resample;
% tstart_per_cell(tstart_per_cell<0) = 0;
% tstart_per_cell = tstart_per_cell/resample;
% 
% 
%% Finish assignments
resample = analysis_out.resample_rate;
tstart = floor(tstart*resample)/resample;
analysis_out.maxsum_tstart = repmat(tstart, 1, size(analysis_out.cellnums,2)); %tstart_per_cell;
if filter_on == true
    analysis_out.maxsum_tstart_basis = 'filtered';
else
    analysis_out.maxsum_tstart_basis = 'raw';
end

analysis_out.version.maxsum_tstart = version;

end