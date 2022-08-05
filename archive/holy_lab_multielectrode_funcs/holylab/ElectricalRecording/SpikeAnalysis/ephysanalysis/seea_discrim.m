function analysis_out = seea_discrim(ephysin, analysis_in)
%SEEA_DISCRIM calculates d' (discriminability) index from ephys or analysis struct
% This function works via 'seea' to calculate d' of a cell for
% all valve openings, over intervals (based on delta_r ONLY @ this time)
% 
% analysis_out.delta_r_discrim:  will contain the same subfields as delta_r
%                                (i.e. s0e1...) each containing a
%                                [n_valves x n_valves] matrix with d'
%                                comparisons between each valve pair
% analysis_out.discrim_version:  The version of this function used to
%                                produce the values in analysis_out.delta_r_discrim
%
% Copyright 2008 Julian P. Meeks (Timothy Holy Laboratory)
% 
% Revision History:
% 2008_05_11- Wrote it (JPM)

%% Check

% Versioning:
disc_version = 1.0;

if isfield(analysis_in,'seea_discrim_version')
    if analysis_in.seea_discrim_version == disc_version
        fprintf('%s seea_discrim is up to date. Version %.1f\n', analysis_in.id, analysis_in.seea_discrim_version);
        analysis_out = analysis_in;
    else
        fprintf('%s seea_discrim is out of date. Updating from Version %.1f to %.1f\n', analysis_in.id, analysis_in.seea_discrim_version, disc_version);
    end
end

analysis_out = analysis_in;

% check for seea_deltar being done
if ~isfield(analysis_out,'delta_r')
    analysis_out = seea_deltar(ephysin, analysis_out);
end

%% Set up seea_discrim matrix
n_valves = size(analysis_out.valvelabels,2);
discrim = zeros(n_valves);                  % this will be a [n_valves x n_valves] matrix
field_names = fieldnames(analysis_out.delta_r);
fieldname_size = size(field_names,1);
n_cells = size(analysis_out.delta_r_mean.(field_names{1}),2);

%% Calculate d' @ each valve, interval, cell

for c_idx = 1:n_cells
    for f_idx = 1:fieldname_size
        for v_idx = 1:n_valves
            thisfield = field_names{f_idx};
            for vv_idx = 1:n_valves
                analysis_out.delta_r_discrim.(field_names{f_idx})(v_idx,vv_idx) =...
                                  dprime(analysis_out.delta_r.(field_names{f_idx}){v_idx,c_idx},...
                                         analysis_out.delta_r.(field_names{f_idx}){vv_idx,c_idx});
            end
        end
    end
end

analysis_out.seea_discrim_version = disc_version;

end