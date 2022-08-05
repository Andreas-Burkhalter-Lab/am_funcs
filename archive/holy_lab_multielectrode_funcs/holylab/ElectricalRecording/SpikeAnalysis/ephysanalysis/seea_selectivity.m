function analysis_out = seea_selectivity(ephysin, analysis_in)
% seea_selectivity calculates gender selectivity on all
% concentrations (dilutions) of male and female urine in a given
% experiment.  If the experiment does not contain both male and female
% urine at equivalent concentrations, no selectivity index will be calculated.
% The input analysis_in will have fields:
%           fem_selectivity_mean
%           fem_selectivity_sem
%           fem_selectivity_ttest
%           fem_selectivity_dilution
%           fem_selectivity_f_valves
%           fem_selectivity_m_valves
% added as fields within the analysis structure.
% NOTE: This function requires analysis.delta_r and/or
% analysis.delta_r_mean.  If fields don't exist the seea_deltar() function
% will be called.
%
% Copyright 2008 Julian P Meeks (Timothy Holy Laboratory)

% Revision History:
% 2008_05_01-2008_05_11: Wrote backbone (JPM)

%% Check data
% Versioning:
sel_version = 1.0; 
if isfield(analysis_in, 'seea_selectivity_version')
    if analysis_in.seea_selectivity_version == sel_version;
        fprintf('%s seea_selectivity is up to date. Version %.1f\n', ephysin(1).basefilename,analysis_in.seea_selectivity_version);
        analysis_out = analysis_in;
        return;
    else
        fprintf('%s seea_selectivity is out of date. Upgrading from version %.1f to version %.1f\n', ephysin(1).basefilename,analysis_in.seea_selectivity_version, sel_version);        
    end
end

analysis_out = analysis_in;

% first check if field already exists.  If so, do nothing special (for now)
if isfield(analysis_out, 'fem_selectivity_mean')
    analysis_out.fem_selectivity_mean = struct;
else
    analysis_out.fem_selectivity_mean = struct;
end
if isfield(analysis_out, 'fem_selectivity_sem')
    analysis_out.fem_selectivity_sem = struct;
else
    analysis_out.fem_selectivity_sem = struct;
end
if isfield(analysis_out, 'fem_selectivity_signif')
    analysis_out.fem_selectivity_signif = struct;
else
    analysis_out.fem_selectivity_signif = struct;
end

% check to make sure delta_r_mean is present or can be calculated
if ~isfield(analysis_in, 'delta_r_mean')&&~isfield(analysis_in, 'delta_r')
    analysis_out = seea_deltar(ephysin, analysis_out);
end

%% find concentrations (dilutions) of male and female urine using valvelabels
male_valves = []; female_valves = [];

% find the valves which have only male or female urine (requires JPM
% standard valve notation MU = male urine, FU = female urine
if ~isfield(analysis_out, 'valvelabels')
    analysis_out.valvelabels = ephysin(1).valvelabels;
end
valvelabels = analysis_out.valvelabels;

for valve_idx = 1:size(valvelabels,2)
    if ~isempty(strfind(valvelabels{valve_idx}, 'MU'))... 
       && isempty(strfind(valvelabels{valve_idx}, '/'))
        if isempty(male_valves)
            male_valves = valve_idx;
        else
            male_valves(end+1) = valve_idx;
        end
    end
    if ~isempty(strfind(valvelabels{valve_idx}, 'FU'))... 
       && isempty(strfind(valvelabels{valve_idx}, '/'))
        if isempty(female_valves)
            female_valves = valve_idx;
        else
            female_valves(end+1) = valve_idx;
        end
    end
end

% parse the dilutions of urine used
% male:
for m_idx = 1:size(male_valves,2)
    coln = strfind(valvelabels{male_valves(m_idx)}, ':');
    spac = strfind(valvelabels{male_valves(m_idx)}, ' ');
    m_dilution(m_idx) = str2num(valvelabels{male_valves(m_idx)}(coln+1:spac-1));
end
for f_idx = 1:size(female_valves,2)
    coln = strfind(valvelabels{female_valves(f_idx)}, ':');
    spac = strfind(valvelabels{female_valves(f_idx)}, ' ');
    f_dilution(f_idx) = str2num(valvelabels{female_valves(f_idx)}(coln+1:spac-1));
end

% choose dilutions which intersect:
[fem_selectivity_dilution,ia, ib] = intersect(m_dilution, f_dilution);

analysis_out.fem_selectivity_f_valves = female_valves;
analysis_out.fem_selectivity_m_valves = male_valves;

% determine fields and n_cells of analysis struct to cycle through:
% first, check for 'analysis_out.defined_subintervals' for which to do the
% analysis:
if isfield(analysis_out, 'defined_subintervals')
    if iscell(analysis_out.defined_subintervals)
        field_names = analysis_out.defined_subintervals;
    elseif analysis_out.defined_subintevals == 'all'
        field_names = fieldnames(analysis_out.delta_r_mean);
        field_names = field_names';
    else
        field_names = {analysis.defined_subintervals};
    end
end
fieldname_size = size(field_names,2);
n_cells = size(analysis_out.delta_r_mean.(field_names{1}),2);
if find(strcmp(analysis_out.valvelabels,'Ringer''s'));
    ringer_valve = find(strcmp(analysis_out.valvelabels,'Ringer''s'));
else
    ringer_valve = [];
end

% calculate selectivity, stdev, sem, ttest, ranksum @ each dilution
for c_idx = 1:n_cells
    for f_idx = 1:fieldname_size
        for d_idx = 1:size(fem_selectivity_dilution,2)
            thisfield = field_names{f_idx};
            this_f_valve = female_valves(ib(d_idx));
            this_m_valve = male_valves(ia(d_idx));
            if ~isempty(ringer_valve)
                f_sample = analysis_out.delta_r.(thisfield){this_f_valve,c_idx};
                m_sample = analysis_out.delta_r.(thisfield){this_m_valve,c_idx};
                r_sample = analysis_out.delta_r.(thisfield){ringer_valve,c_idx};
                if sum(r_sample) == 0;
                    r_sample(1:2:end) = 0.01;
                    r_sample(2:2:end) = -0.01;
                end
            else
                % in the case no Ringer's valve was used, take the firing rate versus the 
                % base firing rate from female valve openings instead (not ideal but hey)
                f_sample = analysis_out.fr.(thisfield){this_f_valve,c_idx};
                m_sample = analysis_out.fr.(thisfield){this_m_valve,c_idx};
                r_sample = analysis_out.base_r{this_f_valve,c_idx};
                if std(r_sample) == 0;
                    r_sample(1:2:end) = mean(r_sample(1:2:end))+0.01;
                    r_sample(2:2:end) = mean(r_sample(2:2:end))-0.01;
                end
                if std(m_sample) == 0;
                    m_sample(1:2:end) = mean(m_sample(1:2:end))+0.01;
                    m_sample(2:2:end) = mean(m_sample(2:2:end))-0.01;
                end
                if std(f_sample) == 0;
                    f_sample(1:2:end) = mean(f_sample(1:2:end))+0.01;
                    f_sample(2:2:end) = mean(f_sample(2:2:end))-0.01;
                end
            end
            analysis_out.fem_selectivity_mean.(thisfield)(d_idx,c_idx) = ...
                calc_fem_selectivity(f_sample,...
                                     m_sample,...
                                     r_sample);
            % for now, compute standard error and t_test p-value using
            % bootstrapping 
          % NOTE: First, must ensure each sample has same number of
          % trials (otherwise bootstrp gets angry)
            size_f = size(f_sample,2);
            size_m = size(m_sample,2);
            size_r = size(r_sample,2);
            min_size = min([size_f size_m size_r]);
            
            f_s_bootstats = bootstrp(10000, @calc_fem_selectivity,...
                rot90(f_sample(1:min_size)),...
                rot90(m_sample(1:min_size)),...
                rot90(r_sample(1:min_size)));
            r_s_bootstats = bootstrp(10000, @calc_fem_selectivity,...
                rot90(r_sample(1:min_size)),...
                flipdim(rot90(r_sample(1:min_size)),1),...
                rot90(r_sample(1:min_size)));
            % calculate selectivity vs. ringers for m/f urine
            if size(m_sample,2)>2 
                cim = bootci(10000, @calc_fem_selectivity,...
                    rot90(m_sample(1:min_size)),...
                    rot90(r_sample(1:min_size)),...
                    rot90(r_sample(1:min_size)));
            else
                cim = 0;
            end
            if size(f_sample,2)>2
                cif = bootci(10000, @calc_fem_selectivity,...
                    rot90(f_sample(1:min_size)),...
                    rot90(r_sample(1:min_size)),...
                    rot90(r_sample(1:min_size)));
            else
                cif = 0;
            end
            analysis_out.fem_selectivity_sem.(thisfield)(d_idx,c_idx) = std(f_s_bootstats);
            if (min(cim)>=0 && max(cim)>=0) || (max(cim)<0 && min(cim)<0)
                passed = 1;
            elseif (min(cif)>=0 && max(cif)>=0) || (max(cif)<0 && min(cif)<0)
                passed = 1;
            else
                passed = 0;
            end
            analysis_out.fem_selectivity_signif.(thisfield)(d_idx,c_idx) = passed;
            
        end
    end
end

analysis_out.fem_selectivity_dilution = fem_selectivity_dilution;
analysis_out.seea_selectivity_version = sel_version;

end