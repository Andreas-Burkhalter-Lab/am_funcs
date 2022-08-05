function analysis_out = seea_lifetime_sparseness(ephysin, analysis_in)
% seea_lifetime_sparseness calculates lifetime sparseness for sulfated
% steroids in a given experiment.  If the experiment does not contain at
% least one significant response to sulfated steroids,
% seea_lifetime_sparseness_version will be set to (current version) and an
% empty matrix returned in analysis.lifetime_sparseness areas
% The input analysis_in will have fields:
%           lifetime_sparseness_version
%           lifetime_sparseness_concentration (in uM)
%           lifetime_sparseness_compounds (in steraloids compound format)
%           lifetime_sparseness_valvenums (indices for matching up)
%           lifetime_sparseness
%           lifetime_sparseness_sem (via bootstrapping)
%           lifetime_sparseness_comp_rank (rank of preference for
%           compounds)
% added as fields within the analysis structure.
% NOTE: This function requires analysis.delta_r and
% analysis.delta_r_mean.  If fields don't exist the seea_deltar() function
% will be called.
%
% Copyright 2008 Julian P Meeks (Timothy Holy Laboratory)

% Revision History:
% 2008_06_11: Wrote backbone (JPM)

%% Check data
% Versioning:
sparse_version = 1.2; 
if isfield(analysis_in, 'seea_lifetime_sparseness_version')
    if analysis_in.seea_lifetime_sparseness_version == sparse_version;
        fprintf('%s seea_lifetime_sparseness is up to date. Version %.1f\n', ephysin(1).basefilename,analysis_in.seea_lifetime_sparseness_version);
        analysis_out = analysis_in;
        return;
    else
        fprintf('%s seea_lifetime_sparseness is out of date. Upgrading from version %.1f to version %.1f\n', ephysin(1).basefilename,analysis_in.seea_lifetime_sparseness_version, sparse_version);        
    end
end

analysis_out = analysis_in;

% first check if field already exists.  If so, do nothing special (for now)
if isfield(analysis_out, 'lifetime_sparseness')
    analysis_out.lifetime_sparseness = struct;
else
    analysis_out.lifetime_sparseness = struct;
end
if isfield(analysis_out, 'lifetime_sparseness_sem')
    analysis_out.lifetime_sparseness_sem = struct;
else
    analysis_out.lifetime_sparseness_sem = struct;
end
if isfield(analysis_out, 'lifetime_sparseness_comp_rank')
    analysis_out.lifetime_sparseness_comp_rank = struct;
else
    analysis_out.lifetime_sparseness_comp_rank = struct;
end
if isfield(analysis_out, 'lifetime_sparseness_compounds')
    analysis_out.lifetime_sparseness_compounds = cell(1);
else
    analysis_out.lifetime_sparseness_compounds = cell(1);
end
if isfield(analysis_out, 'lifetime_sparseness_concentration')
    analysis_out.lifetime_sparseness_concentration = [];
else
    analysis_out.lifetime_sparseness_concentration = [];
end
if isfield(analysis_out, 'lifetime_sparseness_valvenums')
    analysis_out.lifetime_sparseness_valvenums = [];
else
    analysis_out.lifetime_sparseness_valvenums = [];
end

% check to make sure delta_r_mean is present or can be calculated
if ~isfield(analysis_in, 'delta_r_mean')&&~isfield(analysis_in, 'delta_r')
    analysis_out = seea_deltar(ephysin, analysis_out);
end

%% find concentrations of sulfated steroid compounds using valvelabels

% cycle through and find valves which have sulfated steroids
valvelabels = analysis_out.valvelabels;
compounds = cell(1);
valvenums = [];
concentration = [];
for valve_idx = 1:size(valvelabels,2)
    if ~isempty(valvelabels{valve_idx})
        if~isempty(regexp(valvelabels{valve_idx},'[A E P Q]\d\d\d\d'))
            temp = regexp(valvelabels{valve_idx}, '[A E P Q]\d\d\d\d');
            if isempty(compounds{1})
                compounds{1} = valvelabels{valve_idx}(temp:temp+4);
                valvenums(1) = valve_idx;
            else
                compounds{end+1} = valvelabels{valve_idx}(temp:temp+4);
                valvenums(end+1) = valve_idx;
            end
            if (strfind(valvelabels{valve_idx}, 'nM'))
                temp = strfind(valvelabels{valve_idx}, 'nM');
            elseif (strfind(valvelabels{valve_idx}, 'uM'))
                temp = strfind(valvelabels{valve_idx}, 'uM');;
            elseif (strfind(valvelabels{valve_idx}, 'mM'))
                temp = strfind(valvelabels{valve_idx}, 'mM')
            else
                temp = [];
            end
            if ~isempty(temp)
                if isempty(concentration)
                    concentration(1) = str2num(valvelabels{valve_idx}(1:temp-1));
                else
                    concentration(end+1) = str2num(valvelabels{valve_idx}(1:temp-1));
                end
            end
        end
    end
end
analysis_out.lifetime_sparseness_valvenums = valvenums;
analysis_out.lifetime_sparseness_compounds = compounds;
analysis_out.lifetime_sparseness_concentration = concentration;

% determine number of concentrations
concs = unique(analysis_out.lifetime_sparseness_concentration);
nconcs = size(concs,2);
% TO DO - MAKE VALVENUMS INTO SOMETHING SORTED IN CASE DIFFERENT
% CONCENTRATIONS EXIST!!!!

% determine fieldnames based on seea_delta_r
field_names = fieldnames(analysis_out.delta_r);

% determine number of cells based on seea_delta_r
n_cells = size(analysis_out.delta_r_mean.(field_names{1}),2);

% cycle through and calculate lifetime sparseness
% structure will be: 
% analysis_out.lifetime_sparseness...(s%e%fieldname)(conc#, cell#)
for f_idx = 1:size(field_names,1)
    for conc_idx = 1:nconcs   %%% CURRENTLY NO SUPPORT FOR MULTIPLE CONCS>>> SPARSENESS VALS WILL BE THE SAME
        for cell_idx = 1:n_cells
            delta_r_ttest = analysis_out.delta_r_ttest.(field_names{f_idx})(valvenums,cell_idx);
            delta_r_means = analysis_out.delta_r_mean.(field_names{f_idx})(valvenums,cell_idx);
            if min(delta_r_ttest)>=0.1||max(abs(delta_r_means)) < 0.5
                analysis_out.lifetime_sparseness.(field_names{f_idx})(conc_idx, cell_idx) = NaN;
                analysis_out.lifetime_sparseness_sem.(field_names{f_idx})(conc_idx, cell_idx) = NaN;
                analysis_out.lifetime_sparseness_comp_rank.(field_names{f_idx})(:,conc_idx, cell_idx) = {NaN};
                continue;
            end                
            for sam_idx = 1:size(valvenums,2)
                this_dr = analysis_out.delta_r.(field_names{f_idx}){valvenums(sam_idx),cell_idx};
                this_dr_size = size(this_dr);
                if exist('delta_rs', 'var')
                    main_dr = delta_rs;
                    main_dr_size = size(delta_rs);
                else
                    main_dr_size = this_dr_size;
                    main_dr = this_dr;
                end
                if  main_dr_size(2) > this_dr_size(2)
                    pad = zeros(1,main_dr_size(2) - this_dr_size(2));
                    pad(:) = NaN;
                    this_dr = cat(2, this_dr, pad);
                elseif main_dr_size(2) < this_dr_size(2)
                    pad = zeros(main_dr_size(1), this_dr_size(2)-main_dr_size(2));
                    pad(:) = NaN;
                    main_dr = cat(2, main_dr, pad);
                    delta_rs = main_dr;
                end
                delta_rs(sam_idx,:) = this_dr;
                notnan_size(conc_idx, cell_idx) = this_dr_size(2);
            end
            min_dr_size = min(notnan_size(conc_idx,cell_idx));
            analysis_out.lifetime_sparseness.(field_names{f_idx})(conc_idx, cell_idx) = ...
                calc_lifetime_sparseness(delta_r_means);
            bootstats = bootstrp(10000, @calc_lifetime_sparseness, delta_rs(:,1:min_dr_size));
            analysis_out.lifetime_sparseness_sem.(field_names{f_idx})(conc_idx, cell_idx) = ...
                std(bootstats);
            drm_copy = delta_r_means;
            for sam_idx = 1:size(valvenums,2)
                [temp, comp_rank(sam_idx,1)] = max(drm_copy);
                drm_copy(comp_rank) = -10e10;
            end
            analysis_out.lifetime_sparseness_comp_rank.(field_names{f_idx}){:,conc_idx,cell_idx} = comp_rank;
            comp_rank = [];
        end
    end
end

analysis_out.seea_lifetime_sparseness_version = sparse_version;

end