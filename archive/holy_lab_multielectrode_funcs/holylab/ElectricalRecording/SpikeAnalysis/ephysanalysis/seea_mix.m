function analysis_out = seea_mix(ephysin, analysis_in)
%SEEA_MIX calculates differences between mixtures of identical stimuli
% seea_mix calculates the absolute and relative changes in firing rate when
%  urine dilutions are mixed in equal parts.  
% The input analysis_in will have fields:
%           seea_mix_version : Version number
%           mix_dilution (the -fold dilutions used)
%           mix_pref (the preferred stimulus for this interval)
%           mix_supp (mixture "suppression" from "preferred stimulus")
%            mix_supp_mean
%            mix_supp_sem
%            mix_supp_ttest
%           mix_subadd ("subadditivity")
%            etc. as above
% added as fields within the analysis structure.
% NOTE: This function requires analysis.delta_r and/or
% analysis.delta_r_mean.  If fields don't exist the seea_deltar() function
% will be called.
%
% Copyright 2008 Julian P Meeks (Timothy Holy Laboratory)

% Revision History:
% 2008_05_19-2008_05_??: Wrote backbone (JPM)

%% Check data
% Versioning:
mix_version = 1.1; 
if isfield(analysis_in, 'seea_mix_version')
    if analysis_in.seea_mix_version == mix_version;
        fprintf('%s seea_mix is up to date. Version %.1f\n', ephysin(1).basefilename,analysis_in.seea_mix_version);
        analysis_out = analysis_in;
        return;
    else
        fprintf('%s seea_mix is out of date. Upgrading from version %.1f to version %.1f\n', ephysin(1).basefilename,analysis_in.seea_mix_version, mix_version);        
    end
end

analysis_out = analysis_in;

% check to make sure delta_r_mean is present or can be calculated
if ~isfield(analysis_in, 'delta_r_mean')&&~isfield(analysis_in, 'delta_r')
    analysis_out = seea_deltar(ephysin, analysis_out);
end

% check to make sure seea_selectivity has been run (helpful for dilution
% comparisons)
if ~isfield(analysis_out, 'fem_selectivity_mean')
    analysis_out = seea_selectivity(ephysin, analysis_out);
end
%% Determine whether mixtures are present
valvelabels = analysis_out.valvelabels;
for valve_idx = 1:size(valvelabels,2)
    if ~isempty(strfind(valvelabels{valve_idx}, 'M/FU'))
        if ~exist('mix_valves', 'var')
            mix_valves = valve_idx;
        else
            mix_valves(end+1) = valve_idx;
        end
    end
end
if ~exist('mix_valves', 'var')
    fprintf('%s does not appear to have any mixtures to test. SEEA_MIX abandoned.\n', analysis_in.id);
    return;
end 
%% Initialize new fields

% first check if field already exists.  If so, do nothing special (for now)
if isfield(analysis_out, 'mix_dilution')
    analysis_out.mix_dilution = [];
else
    analysis_out.mix_dilution = [];
end
if isfield(analysis_out, 'mix_valves')
    analysis_out.mix_valves = mix_valves;
else
    analysis_out.mix_valves = mix_valves;
end
if isfield(analysis_out, 'mix_pref')
    analysis_out.mix_pref = struct;
else
    analysis_out.mix_pref = struct;
end
if isfield(analysis_out, 'mix_supp_mean')
    analysis_out.mix_supp_mean = struct;
else
    analysis_out.mix_supp_mean = struct;
end
if isfield(analysis_out, 'mix_supp_sem')
    analysis_out.mix_supp_sem = struct;
else
    analysis_out.mix_supp_sem = struct;
end
if isfield(analysis_out, 'mix_subadd_mean')
    analysis_out.mix_subadd_mean = struct;
else
    analysis_out.mix_subadd_mean = struct;
end
if isfield(analysis_out, 'mix_subadd_sem')
    analysis_out.mix_subadd_sem = struct;
else
    analysis_out.mix_subadd_sem = struct;
end
%% Determine mixture dilutions
for mix_idx = 1:size(mix_valves,2)
    coln = strfind(valvelabels{mix_valves(mix_idx)}, ':');
    spac = strfind(valvelabels{mix_valves(mix_idx)}, ' ');
    mix_dilution(mix_idx) = str2num(valvelabels{mix_valves(mix_idx)}(coln+1:spac-1));
end
analysis_out.mix_dilution = mix_dilution;

%% Calculate subadditivity ".mix_subadd..."
% find matching indices to male/female valves for comparison
% Eqn:  (delta_r_male + delta_r_male - delta_r_mix)
%      /(delta_r_male + delta_r_female + fr_ringer's)
for dil_idx = 1:size(mix_dilution,2)
    selectivity_dil_indices(dil_idx) = find(analysis_out.fem_selectivity_dilution == mix_dilution(dil_idx));
end
female_valves = analysis_out.fem_selectivity_f_valves(selectivity_dil_indices);
male_valves = analysis_out.fem_selectivity_m_valves(selectivity_dil_indices);

% find ringer's valve (exclude samples without for now)
if find(strcmp(analysis_out.valvelabels,'Ringer''s'));
    ringer_valve = find(strcmp(analysis_out.valvelabels,'Ringer''s'));
else
    ringer_valve = [];
    fprintf('No Ringer''s control valve found for %s.  SEEA_MIX abandoned.', analysis_out.id);
    return;
end

% Set up for looping
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

% Now, loop through cells, (subinterval) fields, dilutions and calculate
% subadditivity.  
for c_idx = 1:n_cells
    for f_idx = 1:fieldname_size
        for d_idx = 1:size(mix_dilution,2)
            thisfield = field_names{f_idx};
            this_fem_valve = female_valves(d_idx);
            this_mal_valve = male_valves(d_idx);
            this_mix_valve = mix_valves(d_idx);
            
            fem_sample = analysis_out.delta_r.(thisfield){this_fem_valve,c_idx};
            mal_sample = analysis_out.delta_r.(thisfield){this_mal_valve,c_idx};
            mix_sample = analysis_out.delta_r.(thisfield){this_mix_valve,c_idx};
            r_sample = analysis_out.fr.(thisfield){ringer_valve,c_idx};
            fem_sample_mean = analysis_out.delta_r_mean.(thisfield)(this_fem_valve,c_idx);
            mal_sample_mean = analysis_out.delta_r_mean.(thisfield)(this_mal_valve,c_idx);
            mix_sample_mean = analysis_out.delta_r_mean.(thisfield)(this_mix_valve,c_idx);
            r_sample_mean = analysis_out.fr_mean.(thisfield)(ringer_valve,c_idx);
            if sum(r_sample) == 0;
                r_sample(1:2:end) = 0.01;
                r_sample(2:2:end) = -0.01;
            end
            
            % here it is, the equation!
            analysis_out.mix_subadd_mean.(thisfield)(d_idx,c_idx) = ...
                 calc_mix_subadd(fem_sample_mean, mal_sample_mean, mix_sample_mean,r_sample_mean);
            % to find the standard error of the metric: fun with bootstrapping!
            % NOTE: First, must ensure each sample has same number of
            % trials (otherwise bootstrp gets angry)
            size_fem = size(fem_sample,2);
            size_mal = size(mal_sample,2);
            size_mix = size(mix_sample,2);
            size_r = size(r_sample,2);
            min_size = min([size_fem size_mal size_mix size_r]);
            
           mix_bootstats = bootstrp(10000, @calc_mix_subadd,...
               rot90(fem_sample(1:min_size)),...
               rot90(mal_sample(1:min_size)),...
               rot90(mix_sample(1:min_size)),...
               rot90(r_sample(1:min_size)));
           analysis_out.mix_subadd_sem.(thisfield)(d_idx,c_idx) = std(mix_bootstats);
        end
    end
end

%% Calculate mixture suppression ".mix_supp..."
% Eqn:  (delta_r_mix - delta_r_preferred)
%      /(delta_r_preferred + fr_ringer's)

for c_idx = 1:n_cells
    for f_idx = 1:fieldname_size
        for d_idx = 1:size(mix_dilution,2)
            % determine preferred stimulus based on the fem_selectivity
            % NOTE: Will only calculate on intervals where fem_selectivity
            % signif is 1 (true) NaN otherwise
            thisfield = field_names{f_idx};
            this_fem_valve = female_valves(d_idx);
            this_mal_valve = male_valves(d_idx);
            this_mix_valve = mix_valves(d_idx);
            if analysis_out.fem_selectivity_signif.(thisfield)(selectivity_dil_indices(d_idx),c_idx) == 1
              if analysis_out.fem_selectivity_mean.(thisfield)(selectivity_dil_indices(d_idx),c_idx) < 0
                  analysis_out.mix_pref.(thisfield)(d_idx,c_idx) = 'M';
                  this_pref_valve = this_mal_valve;
              elseif analysis_out.fem_selectivity_mean.(thisfield)(selectivity_dil_indices(d_idx),c_idx) > 0
                  analysis_out.mix_pref.(thisfield)(d_idx,c_idx) = 'F';
                  this_pref_valve = this_fem_valve;
              end
            else
                analyis_out.mix_pref.(thisfield)(d_idx,c_idx) = NaN;
                this_pref_valve = [];
            end
            
            if isempty(this_pref_valve)
               analysis_out.mix_supp_mean.(thisfield)(d_idx,c_idx) = NaN;
               analysis_out.mix_supp_sem.(thisfield)(d_idx,c_idx) = NaN;
            else
                pref_sample = analysis_out.delta_r.(thisfield){this_pref_valve,c_idx};
                mix_sample = analysis_out.delta_r.(thisfield){this_mix_valve,c_idx};
                r_sample = analysis_out.fr.(thisfield){ringer_valve,c_idx};
                pref_sample_mean = analysis_out.delta_r_mean.(thisfield)(this_pref_valve,c_idx);
                mix_sample_mean = analysis_out.delta_r_mean.(thisfield)(this_mix_valve,c_idx);
                r_sample_mean = analysis_out.fr_mean.(thisfield)(ringer_valve,c_idx);
                % now calculate mix_supp_mean!!
                analysis_out.mix_supp_mean.(thisfield)(d_idx,c_idx) = ...
                    calc_mix_supp(pref_sample_mean, mix_sample_mean, r_sample_mean);
                mix_bootstats = bootstrp(10000, @calc_mix_supp,...
                    rot90(pref_sample(1:min_size)),...
                    rot90(mix_sample(1:min_size)),...
                    rot90(r_sample(1:min_size)));
                analysis_out.mix_supp_sem.(thisfield)(d_idx,c_idx) = std(mix_bootstats);
            end         
        end
    end
end

analysis_out.seea_mix_version = mix_version;

end