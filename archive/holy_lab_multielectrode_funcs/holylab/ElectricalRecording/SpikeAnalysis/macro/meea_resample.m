function analysis_out = meea_resample(ephysin, analysis_in)
% meea_resample adds a resampled firing rate at the desired rate(s)
% Syntax: analysis_out = meea_resample(ephysin, analysis_in)
%
% Inputs:  ephysin: multi electrode ephys structure
%          analysis_in: analysis structure with optional subfields:
%            trange: scalar or vector of desired region(s) to resample (must be within ephysin trange to function)
%            resample_rate: scalar, desired resampling rate Hz
%            resample_filter: 1 x 2 matrix: [b a] 1-D filter (default = cheby1(4,0.5,resample_rate/5)) 
% Outputs: analysis_out structure variable will be a copy of analysis in with the following fields added
%            .trange (if not present, default will be trange of ephysin)
%            .resample_rate: if not present, default will be 10 Hz)
%            .resample_filter: if not present, default will be cheby1(4, 0.5, 0.2) = 1 Hz LPF @ 10 Hz resample
%            .resampled_raw: {{valve #}(trial, samples, cell)}
%            .resampled_filter: {{valve #}(trial, samples, cell)}
%            .verstion.resample: current version number
% See also multi_electrode_ephys_analyze, ephys, cheby1

% Copyright 2009 Julian P. Meeks (Timothy Holy Laboratory)
%
% Version History:
% 2009_01-2009_05: Wrote it (JPM)
% 2009_08_13: checked and fixed help info (JPM)
%% Version control
version = 1.02;
if isfield(analysis_in, 'version')
    if isfield(analysis_in.version, 'resample')
        if analysis_in.version.resample == version
            fprintf('meea_resample is up to date.\n');
            analysis_out = analysis_in;
            return;
        else
            fprintf('meea_resample is out of date. Updating from v. %2d to v. %2d.\n', analysis_in.version.resample, version);
        end
    end
end

%% initialize default variables
analysis_out = analysis_in;
analysis_out = default(analysis_out, 'resample_rate', 10);
[b, a] = cheby1(4, 0.5, 0.2);
analysis_out = default(analysis_out, 'resample_filter', [b;a]);

%% Fxn main
rate = analysis_out.resample_rate;
ncells = size(ephysin(1).cellnums,2);
wholematrix = [];
% if filter supplied via options.filter_ba, set here
filter_on = false;
if isfield(analysis_out, 'resample_filter')
    if ~isempty(analysis_out.resample_filter)
        b = analysis_out.resample_filter(1,:);
        a = analysis_out.resample_filter(2,:);
        filter_on = true;
    end
end

% index valve-to-interval
utags = analysis_out.valvelabels;
ntags = size(utags,2);
valve2interval = cell(1);
% cycle per valve and strip out intervals
for vi = 1:ntags
    count = 1;
    for di = 1:size(ephysin,2)
        if ~isempty(strmatch(utags{vi},ephysin(di).tag, 'exact'))
            valve2interval{vi}(count) = di;
            count = count +1;
        end
    end
end
% go through each interval and resample
trange = [ephysin(1).toffset diff(ephysin(1).scanrange/ephysin(1).scanrate)+ephysin(1).toffset];
for vi = 1:ntags
    nseconds = round(diff(trange));
    hz_desired = rate;
    this_matrix = zeros([size(valve2interval{vi},2) hz_desired*nseconds-1 ncells]);
    for di = 1:size(valve2interval{vi},2)
        for ci = 1:ncells
            thisinterval = valve2interval{vi}(di);
            thesespikes = ephysin(thisinterval).celltimes{ci};
            desired_scan_bin = (ephysin(thisinterval).scanrate/hz_desired);
            cnums(ci)=ephysin(thisinterval).cellnums(ci);
            for ri = 1:hz_desired*nseconds-1
                thisstart = ephysin(thisinterval).scanrange(1)+(ri-1)*desired_scan_bin;
                thisend = thisstart+desired_scan_bin;
                thisrange = [thisstart thisend];
                nspikes_in_bin = nnz(isbetween(thesespikes, thisrange));
                if isempty(nspikes_in_bin)
                    this_matrix(di,ri,ci) = 0;
                else
                    this_matrix(di,ri,ci) = nspikes_in_bin*hz_desired;
                end
            end
        end
    end
    if filter_on
        filtmatrix{vi} = filter(b,a,this_matrix,[],2);
    end
    wholematrix{vi} = this_matrix(:,:,:);
end

%% Finish assignments
analysis_out.resampled_raw = wholematrix;
analysis_out.resampled_filter = filtmatrix;
analysis_out.trange = trange;
analysis_out.cellnums = cnums;

analysis_out.version.resample = version;


end