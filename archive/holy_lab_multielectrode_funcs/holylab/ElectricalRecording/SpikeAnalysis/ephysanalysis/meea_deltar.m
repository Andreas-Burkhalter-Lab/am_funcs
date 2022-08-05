function analysis_out = meea_deltar(ephysin, analysis_in)
% meea_deltar calculates delta_r values over the designated time range(s)
% Syntax: analysis_out = meea_deltar(ephysin, analysis_in)
%
% ** NOTE: THIS VERSION OF ANALYSIS IS NOT PREFERRED FOR VNO DATA **
%
% Inputs:  ephysin: multi electrode ephys structure
%          analysis_in: analysis structure with required subfields:
%            
% Outputs: analysis_out structure variable will be a copy of analysis in with the following fields added
%            .seea_deltar_version: 1.0, etc... for future use
%            .delta_r_source ('celltimes' or 'sniptimes'; celltimes default)
%            .delta_r_range  (cell array of [start end] vectors in seconds)
%            .delta_r        (for each range, for each trial)
%            .delta_r_mean   (for each range)
%            .delta_r_ttest  (for each range)
%            .delta_r_ranksum (for each range)
%            .fr ...         (raw firing rates and derivatives as with delta_r)
%            .base_r         (average firing rate for 5 seconds preceding
%                             stimulus and derivatives)
%
% See also MULTI_ELECTRODE_EPHYS_ANALYZE SEEA_DELTAR

% % Copyright 2009 Julian P. Meeks (Timothy Holy Laboratory)
% 
% Version History:
% 2009_01-2009_05: Wrote it (JPM)
% 2009_08_13: checked and fixed help info (JPM)

%% check for existing delta_r info, prompt for overwrite vs. append

% Version control

% CURRENT VERSION:
version = 1.2;  % 2008_06_11

if isfield(analysis_in, 'seea_deltar_version')
    if analysis_in.seea_deltar_version == version
        analysis_out = analysis_in;
        fprintf('%s seea_deltar is up to date. Version %.1f\n', ephysin(1).basefilename,analysis_in.seea_deltar_version);
        return;
    else
        fprintf('%s deltar is out of date. Upgrading from version %.1f to version %.1f\n', ephysin(1).basefilename,analysis_in.seea_deltar_version, version);        
    end
end
%

if isfield(analysis_in,'delta_r')
    response = questdlg(sprintf('analysis.delta_r and associated fields exist for %s.  Overwrite (o) or append (a)?:', ephysin(1).basefilename),...
        'Overwrite analysis.delta_r?', '(o)verwrite', '(a)ppend', 'Cancel', '(a)ppend');
    switch response
        case '(o)verwrite'
            analysis_out = analysis_in;
            if isfield(analysis_out,'base_r'); analysis_out = rmfield(analysis_out,'base_r');end;
            if isfield(analysis_out,'base_r_mean'); analysis_out = rmfield(analysis_out,'base_r_mean');end;
            if isfield(analysis_out,'base_r_stdev'); analysis_out = rmfield(analysis_out,'base_r_stdev');end;
            if isfield(analysis_out,'base_r_sem'); analysis_out = rmfield(analysis_out,'base_r_sem');end;
            if isfield(analysis_out,'fr'); analysis_out = rmfield(analysis_out,'fr');end;
            if isfield(analysis_out,'fr'); analysis_out = rmfield(analysis_out,'fr');end;
            if isfield(analysis_out,'fr'); analysis_out = rmfield(analysis_out,'fr');end;
            if isfield(analysis_out,'fr'); analysis_out = rmfield(analysis_out,'fr');end;
            if isfield(analysis_out,'delta_r'); analysis_out = rmfield(analysis_out,'delta_r');end;
            if isfield(analysis_out,'delta_r'); analysis_out = rmfield(analysis_out,'delta_r');end;
            if isfield(analysis_out,'delta_r'); analysis_out = rmfield(analysis_out,'delta_r');end;
            if isfield(analysis_out,'delta_r'); analysis_out = rmfield(analysis_out,'delta_r');end;
        case '(a)ppend'
            analysis_out = analysis_in;
        case 'Cancel'
            analysis_out = analysis_in; return;
    end
else
    analysis_out = analysis_in;
end

%% calculate base_r each epoch in ephysin

for ephys_idx = 1:size(ephysin,2)
    
    % Step 1: determine which valve is open during this interval
    % Important: if duplicate valves have been used, reassign valve numbers
    % to combine like valves
    uniquevalvelabels = unique(ephysin(ephys_idx).valvelabels);
    valvenum = find(strcmp(uniquevalvelabels,ephysin(ephys_idx).tag));
    
    % Step 2: determine whether we are calculating based on celltimes or sniptimes for this interval:
    if isfield(ephysin(ephys_idx), 'celltimes')
        choice = 'celltimes';
    elseif isfield(ephysin(ephys_idx), 'sniptimes')
        choice = 'sniptimes';
    else
        warnmsg = sprintf('%s has not been snippeted.  Delta_r analysis not performed.',ephysin(ephys_idx).basefilename);
        warning(warnmsg);
    end
    % Assign delta_r_source{valvenum}{repeat}
    if ~isfield(analysis_out, 'delta_r_source')
        analysis_out.delta_r_source = cell(size(ephysin(1).valvelabels,2),1);
    end
    if isempty(analysis_out.delta_r_source{valvenum(1)})
        analysis_out.delta_r_source{valvenum(1)}{1} = choice;
    else
        analysis_out.delta_r_source{valvenum(1)}{end+1} = choice;
    end
    
    % Step 3: Initialize variables to be calculated
    
    valve_size = size(unique(ephysin(1).valvelabels),2);
    n_cells = size(ephysin(ephys_idx).(choice),2);
    % 3a: base_r and derivatives for spont. rate in 5 seconds preceding stim
    if ~isfield(analysis_out,'base_r')
        analysis_out.base_r = cell(valve_size, n_cells);
    end
    if ~isfield(analysis_out,'base_r_mean')
        analysis_out.base_r_mean = zeros(valve_size, n_cells);
    end
    if ~isfield(analysis_out,'base_r_stdev')
        analysis_out.base_r_stdev = zeros(valve_size, n_cells);
    end
    if ~isfield(analysis_out, 'base_r_sem')
        analysis_out.base_r_sem = zeros(valve_size, n_cells);
    end
    if ~isfield(analysis_out, 'base_r_ttest')
        analysis_out.base_r_ttest = zeros(valve_size, n_cells);
    end
    if ~isfield(analysis_out, 'base_r_ranksum')
        analysis_out.base_r_ranksum = zeros(valve_size, n_cells);
    end
    % 3b: fr and derivatives for firing rate in varying intervals
    if ~isfield(analysis_out,'fr')
        analysis_out.fr = struct;
    end
    if ~isfield(analysis_out,'fr_mean')
        analysis_out.fr_mean = struct;
    end
    if ~isfield(analysis_out,'fr_stdev')
        analysis_out.fr_stdev = struct;
    end
    if ~isfield(analysis_out, 'fr_sem')
        analysis_out.fr_sem = struct;
    end
    if ~isfield(analysis_out, 'fr_ttest')
        analysis_out.fr_ttest = struct;
    end
    if ~isfield(analysis_out, 'fr_ranksum')
        analysis_out.fr_ranksum = struct;
    end    
    % 3c: delta_r and derivatives for change in firing rate in varying intervals
    if ~isfield(analysis_out,'delta_r')
        analysis_out.delta_r = struct;
    end
    if ~isfield(analysis_out,'delta_r_mean')
        analysis_out.delta_r_mean = struct;
    end
    if ~isfield(analysis_out,'delta_r_stdev')
        analysis_out.delta_r_stdev = struct;
    end
    if ~isfield(analysis_out, 'delta_r_sem')
        analysis_out.delta_r_sem = struct;
    end    
    if ~isfield(analysis_out, 'delta_r_ttest')
        analysis_out.delta_r_ttest = struct;
    end
    if ~isfield(analysis_out, 'delta_r_ranksum')
        analysis_out.delta_r_ranksum = struct;
    end
    
    % Step 4: Step through ephysin
    for cell_idx = 1:size(ephysin(ephys_idx).(choice),2)
        % Check to make sure there are times associated with this interval,
        % if not, set all per-interval/repeat data to 0
        if isempty(ephysin(ephys_idx).(choice){cell_idx})
            ephysin(ephys_idx).(choice){cell_idx} = 0;
        end
        
        if ~isempty(ephysin(ephys_idx).(choice){cell_idx})
        
            % Step 4a: Calculate base_r for this interval
            baserate= size(intersect(ephysin(ephys_idx).(choice){cell_idx}, ...
                [ephysin(ephys_idx).scanrange(1)+(-5-ephysin(ephys_idx).toffset)*ephysin(ephys_idx).scanrate: ...
                ephysin(ephys_idx).scanrange(1)-(ephysin(ephys_idx).toffset*ephysin(ephys_idx).scanrate)]),2)/5;
            if isempty(baserate)
                baserate = 0;
            end
            if isempty(analysis_out.base_r{valvenum(1),cell_idx})
                analysis_out.base_r{valvenum(1),cell_idx}=baserate;
            else
                analysis_out.base_r{valvenum(1),cell_idx}(end+1)=baserate;
            end
        end
    end
end
%% calculate delta_r, base_r, and fr data for each interval in ephysin

for ephys_idx = 1:size(ephysin,2)
    
    % Step 1: determine which valve is open during this interval
    % Important: if duplicate valves have been used, reassign valve numbers
    % to combine like valves
    uniquevalvelabels = unique(ephysin(ephys_idx).valvelabels);
    valvenum = find(strcmp(uniquevalvelabels,ephysin(ephys_idx).tag));
    
    % Step 2: determine whether we are calculating based on celltimes or sniptimes for this interval:
    if isfield(ephysin(ephys_idx), 'celltimes')
        choice = 'celltimes';
    elseif isfield(ephysin(ephys_idx), 'sniptimes')
        choice = 'sniptimes';
    else
        warnmsg = sprintf('%s has not been snippeted.  Delta_r analysis not performed.',ephysin(ephys_idx).basefilename);
        warning(warnmsg);
    end
    % Assign delta_r_source{valvenum}{repeat}
    if ~isfield(analysis_out, 'delta_r_source')
        analysis_out.delta_r_source = cell(size(ephysin(1).valvelabels,2),1);
    end
    if isempty(analysis_out.delta_r_source{valvenum(1)})
        analysis_out.delta_r_source{valvenum(1)}{1} = choice;
    else
        analysis_out.delta_r_source{valvenum(1)}{end+1} = choice;
    end
    
    % Step 3: Initialize variables to be calculated
    
    valve_size = size(unique(ephysin(1).valvelabels),2);
    n_cells = size(ephysin(ephys_idx).(choice),2);
    % 3a: base_r and derivatives for spont. rate in 5 seconds preceding stim
    if ~isfield(analysis_out,'base_r')
        analysis_out.base_r = cell(valve_size, n_cells);
    end
    if ~isfield(analysis_out,'base_r_mean')
        analysis_out.base_r_mean = zeros(valve_size, n_cells);
    end
    if ~isfield(analysis_out,'base_r_stdev')
        analysis_out.base_r_stdev = zeros(valve_size, n_cells);
    end
    if ~isfield(analysis_out, 'base_r_sem')
        analysis_out.base_r_sem = zeros(valve_size, n_cells);
    end
    if ~isfield(analysis_out, 'base_r_ttest')
        analysis_out.base_r_ttest = zeros(valve_size, n_cells);
    end
    if ~isfield(analysis_out, 'base_r_ranksum')
        analysis_out.base_r_ranksum = zeros(valve_size, n_cells);
    end
    % 3b: fr and derivatives for firing rate in varying intervals
    if ~isfield(analysis_out,'fr')
        analysis_out.fr = struct;
    end
    if ~isfield(analysis_out,'fr_mean')
        analysis_out.fr_mean = struct;
    end
    if ~isfield(analysis_out,'fr_stdev')
        analysis_out.fr_stdev = struct;
    end
    if ~isfield(analysis_out, 'fr_sem')
        analysis_out.fr_sem = struct;
    end
    if ~isfield(analysis_out, 'fr_ttest')
        analysis_out.fr_ttest = struct;
    end
    if ~isfield(analysis_out, 'fr_ranksum')
        analysis_out.fr_ranksum = struct;
    end    
    % 3c: delta_r and derivatives for change in firing rate in varying intervals
    if ~isfield(analysis_out,'delta_r')
        analysis_out.delta_r = struct;
    end
    if ~isfield(analysis_out,'delta_r_mean')
        analysis_out.delta_r_mean = struct;
    end
    if ~isfield(analysis_out,'delta_r_stdev')
        analysis_out.delta_r_stdev = struct;
    end
    if ~isfield(analysis_out, 'delta_r_sem')
        analysis_out.delta_r_sem = struct;
    end    
    if ~isfield(analysis_out, 'delta_r_ttest')
        analysis_out.delta_r_ttest = struct;
    end
    if ~isfield(analysis_out, 'delta_r_ranksum')
        analysis_out.delta_r_ranksum = struct;
    end
    
    % Step 4: Step through ephysin
    for cell_idx = 1:size(ephysin(ephys_idx).(choice),2)
        % Check to make sure there are times associated with this interval,
        % if not, set all per-interval/repeat data to 0
        if isempty(ephysin(ephys_idx).(choice){cell_idx})
            ephysin(ephys_idx).(choice){cell_idx} = 0;
        end
        
        if ~isempty(ephysin(ephys_idx).(choice){cell_idx})   
            
            % Step 4b parse 'analysis.defined_subintervals' to choose only
            % the subintervals desired by the user for calculation
            if isfield(analysis_out, 'defined_subintervals')
                if ~iscell(analysis_out.defined_subintervals)
                    if analysis_out.defined_subintervals == 'all'
                        analysis_out.defined_subintervals = ...
                            {'s0e1', 's1e2', 's2e3', 's3e4', 's4e5', 's5e6', 's6e7', 's7e8', 's8e9', 's9e10',...
                            's10e11', 's11e12', 's12e13', 's13e14', 's14e15', 's15e16', 's16e17', 's17e18',...
                            's18e19', 's19e20',...
                            's1e11', 's1e16', 's1e20'};
                        [s_vals e_vals] = strip_se(analysis_out.defined_subintervals);
                    else
                        [s_vals e_vals] = strip_se(analysis_out.defined_subintervals);
                    end
                else
                    [s_vals e_vals] = strip_se(analysis_out.defined_subintervals);
                end
            end

            % use baserate for cell overall, not just for this trial
            temp = 0;
            for vidx = 1:size(analysis_out.valvelabels,2)
                temp = temp+mean(analysis_out.base_r{vidx,cell_idx});
            end
            baserate = temp/size(analysis_out.valvelabels,2);
            analysis_out.mean_baserate = baserate;
            
            % Step 4c: Calculate fr and delta_r for this epoch
            for se_idx = 1:size(s_vals,2)
                if e_vals(se_idx) > s_vals(se_idx)
                    fldname = analysis_out.defined_subintervals{se_idx};
                    thesecells = ephysin(ephys_idx).(choice){cell_idx};
                    thisstart = ephysin(ephys_idx).scanrange(1)+((-ephysin(ephys_idx).toffset)+s_vals(se_idx))*ephysin(ephys_idx).scanrate;
                    thisend = thisstart+(s_vals(se_idx)+e_vals(se_idx))*ephysin(ephys_idx).scanrate;
                    % sr
                    this_sr = ...
                        size(intersect(thesecells, ...
                        [thisstart:thisend]),2)...
                        /(e_vals(se_idx)-s_vals(se_idx));
                    if isempty(this_sr)
                        this_sr = 0;
                    end
                    if ~isfield(analysis_out.fr, fldname)
                        analysis_out.fr.(fldname)=cell(valve_size, n_cells);
                        analysis_out.fr.(fldname){valvenum(1),cell_idx} = this_sr;
                    else
                        analysis_out.fr.(fldname){valvenum(1),cell_idx}(end+1) = this_sr;
                    end
                    % delta_r
                    this_delta_r = ...
                        size(intersect(thesecells, ...
                        [thisstart:thisend]),2)...
                        /(e_vals(se_idx)-s_vals(se_idx))...
                        -baserate;
                    if isempty(this_delta_r)
                        this_delta_r = 0;
                    end
                    if ~isfield(analysis_out.delta_r, fldname)
                        analysis_out.delta_r.(fldname)=cell(valve_size, n_cells);
                        analysis_out.delta_r.(fldname){valvenum(1),cell_idx} = this_delta_r;
                    else
                        analysis_out.delta_r.(fldname){valvenum(1),cell_idx}(end+1) = this_delta_r;
                    end
                end
            end
        end
    end
end
% Step 5: calculate averages, stdevs, sems, ttests, ranksum tests for each valvenum

% for statistical tests, must identify Ringer's valves
if find(strcmp(uniquevalvelabels,'Ringer''s'));
    ringer_valve = find(strcmp(uniquevalvelabels,'Ringer''s'));
else
    ringer_valve = [];
end

% step through each analysis set
for v = 1:valve_size
    for c = 1:n_cells
        if ~isempty(strfind(uniquevalvelabels{v},'EMPTY')) || ~isempty(strfind(uniquevalvelabels{v},'null'))
            continue;
        end
          
        % for base_r
        analysis_out.base_r_mean(v,c)=mean(analysis_out.base_r{v,c});
        analysis_out.base_r_stdev(v,c)=std(analysis_out.base_r{v,c});
        analysis_out.base_r_sem(v,c)=std(analysis_out.base_r{v,c})/sqrt(size(analysis_out.base_r{v,c},2));
        if ~isempty(ringer_valve)
            if size(analysis_out.base_r{v,c},2)>1 && size(analysis_out.base_r{ringer_valve,c},2)>1
                if std(analysis_out.base_r{v,c}) == 0     
                    analysis_out.base_r{v,c} = mean(analysis_out.base_r{v,c})+abs(randn(size(analysis_out.base_r{v,c})))*0.2;
                end
                if std(analysis_out.base_r{ringer_valve,c}) == 0
                    analysis_out.base_r{ringer_valve,c} = mean(analysis_out.base_r{ringer_valve,c})+abs(randn(size(analysis_out.base_r{ringer_valve,c})))*0.2;
                end
                [temp analysis_out.base_r_ttest(v,c)] = ttest2(analysis_out.base_r{v,c},analysis_out.base_r{ringer_valve,c},0.05,'both');
                analysis_out.base_r_ranksum(v,c) = ranksum(analysis_out.base_r{v,c},analysis_out.base_r{ringer_valve,c});
            else
                analysis_out.base_r_ttest(v,c) = NaN;
                analysis_out.base_r_ranksum(v,c) = NaN;
            end
        else
            if isfield(analysis_out, 'base_r_ttest')
                analysis_out = rmfield(analysis_out,'base_r_ttest');
            end
            if isfield(analysis_out,'base_r_ranksum')
                analysis_out = rmfield(analysis_out,'base_r_ranksum');
            end
        end
        temp = [];  % throwaway variable
        % for fr, delta_r, need to step through the individual s%de%d intervals
        for se_idx = 1:size(s_vals,2)
            if e_vals(se_idx) > s_vals(se_idx)
                se_string = sprintf('s%de%d',s_vals(se_idx),e_vals(se_idx));
                analysis_out.fr_mean.(se_string)(v,c)=mean(analysis_out.fr.(se_string){v,c});
                analysis_out.fr_stdev.(se_string)(v,c)=std(analysis_out.fr.(se_string){v,c});
                analysis_out.fr_sem.(se_string)(v,c)=std(analysis_out.fr.(se_string){v,c})/sqrt(size(analysis_out.fr.(se_string){v,c},2));
                if ~isempty(ringer_valve)
                    if size(analysis_out.fr.(se_string){v,c},2)>1 && size(analysis_out.fr.(se_string){ringer_valve,c},2)>1
                        if std(analysis_out.fr.(se_string){v,c}) == 0
                            analysis_out.fr.(se_string){v,c} = abs(mean(analysis_out.fr.(se_string){v,c})+randn(size(analysis_out.fr.(se_string){v,c})))*0.2;
                        end
                        [temp analysis_out.fr_ttest.(se_string)(v,c)] = ttest2(analysis_out.fr.(se_string){v,c},analysis_out.fr.(se_string){ringer_valve,c},0.05,'both');
                        analysis_out.fr_ranksum.(se_string)(v,c) = 2*ranksum(analysis_out.fr.(se_string){v,c},analysis_out.fr.(se_string){ringer_valve,c});
                    else
                        analysis_out.fr_ttest.(se_string)(v,c) = NaN;
                        analysis_out.fr_ranksum.(se_string)(v,c) = NaN;
                    end
                else
                    analysis_out.fr_ttest = NaN;
                    analysis_out.fr_ranksum = NaN;
                end
                analysis_out.delta_r_mean.(se_string)(v,c)=mean(analysis_out.delta_r.(se_string){v,c});
                analysis_out.delta_r_stdev.(se_string)(v,c)=std(analysis_out.delta_r.(se_string){v,c});
                analysis_out.delta_r_sem.(se_string)(v,c)=std(analysis_out.delta_r.(se_string){v,c})/sqrt(size(analysis_out.delta_r.(se_string){v,c},2));
                if ~isempty(ringer_valve)
                    if size(analysis_out.delta_r.(se_string){v,c},2)>1 && size(analysis_out.delta_r.(se_string){ringer_valve,c},2)>1
                        if std(analysis_out.delta_r.(se_string){v,c}) == 0
                            analysis_out.delta_r.(se_string){v,c} = mean(analysis_out.delta_r.(se_string){v,c})+randn(size(analysis_out.delta_r.(se_string){v,c}))*0.2;
                        end
                        [temp analysis_out.delta_r_ttest.(se_string)(v,c)] = ttest2(analysis_out.delta_r.(se_string){v,c},analysis_out.delta_r.(se_string){ringer_valve,c},0.05,'both');
                        analysis_out.delta_r_ranksum.(se_string)(v,c) = 2*ranksum(analysis_out.delta_r.(se_string){v,c},analysis_out.delta_r.(se_string){ringer_valve,c});
                    else
                        analysis_out.delta_r_ttest.(se_string)(v,c) = NaN;
                        analysis_out.delta_r_ranksum.(se_string)(v,c) = NaN;
                    end
           
                else
                    if isfield(analysis_out, 'delta_r_ttest')
                        analysis_out = rmfield(analysis_out,'delta_r_ttest');
                    end
                    if isfield(analysis_out,'delta_r_ranksum')
                        analysis_out = rmfield(analysis_out,'delta_r_ranksum');
                    end
                end
            end
        end
    end
end
% If we've completed this batch without error, stamp the version
analysis_out.seea_deltar_version = version;

end