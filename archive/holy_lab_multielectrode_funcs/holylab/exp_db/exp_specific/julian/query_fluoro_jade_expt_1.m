function matched = query_fluoro_jade_expt_1(entry)
%query_fluoro_jade_expt_1: uses the search_xdb function to sift through
%the exp_db and return the xdb_ids which experiments have been analyzed and
%are ready for plotting
%
% see also: CELLCOUNT_FLAT

% 2007 Julian P. Meeks

% The necessary fields for plotting/analysis are:
%   entry.tags includes 'fluoro_jade'
%   entry.fluoro_jade field exists
%   entry.hist_analysis:
%           -> region_names
%           -> cells_in_regions
%           -> region_areas
%
% In future, will add check that experiment is
%   within the desired dates for this particular
%   experiment
%
% The optional fields
%   entry.hist_analysis
%           -> um_per_pixel
%------------------------------------------------

matched = 0;

if isempty(strmatch('julian', entry.investigator));
    return;
end

if ~isfield(entry, 'tags')
    return;
elseif isempty(strfind(entry.tags, 'fluoro_jade'))
    return;
end 
if ~isfield(entry, 'fluoro_jade'); return; end;

if ~isfield(entry, 'hist_analysis')
    return;
elseif ~isfield(entry.hist_analysis, 'cells_in_regions')...
        || ~isfield(entry.hist_analysis, 'region_names')...
        || ~isfield(entry.hist_analysis, 'region_areas')
    return;
else
    matched = 1;
end

end
