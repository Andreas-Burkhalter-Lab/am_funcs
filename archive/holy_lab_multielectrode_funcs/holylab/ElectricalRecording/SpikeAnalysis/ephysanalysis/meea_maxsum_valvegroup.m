function analysis_out = meea_maxsum_valvegroup(ephysin, analysis_in)
% meea_maxsum_tstart parses valvelabels to group compounds for maxsum_monotonic
%
% syntax: analysis_out = meea_maxsum_valvegroup(ephysin, analysis_in)
%
% inputs:  ephysin: multi electrode ephys structure
%          analysis_in: analysis structure with optional subfields:
%            .valvelabels
% outputs: analysis_out structure variable will be a copy of analysis in with the following fields added
%            .maxsum_valvegroups {ngroups}(n_concentrations in increasing order) - indexes of valveids needed
%            .maxsum_valveids {ngroups} - names of compounds used in .maxsum_valvegroups
% See also multi_electrode_ephys_analyze, ephys, 

% Copyright 2009 Julian P. Meeks (Timothy Holy Laboratory)
% 
%% Version control
version = 1.0;
if isfield(analysis_in, 'version')
    if isfield(analysis_in.version, 'maxsum_valvegroup')
        if analysis_in.version.maxsum_valvegroup == version
            fprintf('meea_maxsum_valvegroup is up to date.\n');
            analysis_out = analysis_in;
            return;
        end
    end
end

%% Fxn main
n_valves = size(analysis_in.valvelabels,2);
remaining_valves = analysis_in.valvelabels;

count = 1;
while ~isempty(remaining_valves);
   thisvalve = remaining_valves{1};
   if strmatch(thisvalve, 'Ringer''s', 'exact')
       valvegroups{count}(1) = strmatch(thisvalve, analysis_in.valvelabels, 'exact');
       valveids{count} = thisvalve;
       count = count+1;
       remaining_valves(1) = [];
   elseif strmatch('1:100', thisvalve)
       valvegroups{count}(1) = strmatch(thisvalve, analysis_in.valvelabels, 'exact');
       valveids{count} = thisvalve;
       count = count+1;
       remaining_valves(1) = [];
   elseif strmatch('ringers', thisvalve)
       valvegroups{count}(1) = strmatch(thisvalve, analysis_in.valvelabels, 'exact');
       valveids{count} = thisvalve;
       count = count+1;
       remaining_valves(1) = [];
   elseif strmatch('vehicle', thisvalve)
       valvegroups{count}(1) = strmatch(thisvalve, analysis_in.valvelabels, 'exact');
       valveids{count} = thisvalve;
       count = count+1;
       remaining_valves(1) = [];   
   elseif strmatch('high', thisvalve)
       valvegroups{count}(1) = strmatch(thisvalve, analysis_in.valvelabels, 'exact');
       valveids{count} = thisvalve;
       count = count+1;
       remaining_valves(1) = [];   
   elseif ~isempty(regexp(thisvalve, '[A E P Q]\d\d\d\d','once'))
       thisid = regexp(thisvalve, '[A E P Q]\d\d\d\d', 'match');
       valveids(count) = thisid;
       same_ids = strmatch(thisid{1}, analysis_in.valvelabels);
       for same_idx = 1:size(same_ids,1)
           thisvalve = analysis_in.valvelabels{same_ids(same_idx)};
           position = 7;
           c = 1;
           tempstr = char;
           while ~isempty(str2num(thisvalve(position))) || ~isempty(strmatch('.', thisvalve(position)))
               tempstr(c) = thisvalve(position);
               c = c+1;
               position = position+1;
           end
           tempnum(same_idx,1) = str2num(tempstr);
           tempnum(same_idx,2) = same_ids(same_idx);           
       end
       [y i] = sort(tempnum(:,1));
       valvegroups{count} = tempnum(i,2)';
       remaining_valves(1:size(same_ids,1)) = [];
       count = count+1;
   end
end

%% Finish assignments
analysis_out = analysis_in;
analysis_out.maxsum_valvegroups = valvegroups;
analysis_out.maxsum_valveids = valveids;

analysis_out.version.maxsum_valvegroup = version;

end