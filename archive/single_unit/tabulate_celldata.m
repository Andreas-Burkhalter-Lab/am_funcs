% create a table with a row containing all data relevant to a given
% cell; run this script after generate_filelist and getTempoTrialData
%%% 9/17/16 AM edited on msi

%% parameters
clear_workspace = 1; % if true, clear all vars except celldata and files_with_dirs after running
%%% vars_to_show lists parameters of interest to be presented in summary
%%% tables
% first column in vars_to_show is the protocol, second is parameter names
% within protocol summary files
vars_to_show = {
%     'RF_MAPPING',{};... % rf mapping with drifting dots, not gratings... not yet included because analysis function doesn't create output results file     
%     'GRATING_RF_MAP',{};... 
    'SIZE_TUNING',... %%% dot or dot field size, not grating size
        {'OptSiz'};...
%     'GRATING_SIZE',{};... 
    'GRATING_SPATIAL_FREQ',...
        {'AnovaP';'MaxLoc';'SFDI'};...
    'GRATING_TEMPORAL_FREQ',...
        {'AnovaP';'MaxLoc';'TFDI'};...
    'SPEED_TUNING',... % drifting dot field speed
        {'AnovaP';'MaxLoc',;'SDI'};...
    'MOTION_COHERENCE',... % drifting dot field motion coherence
        {'PRaw';'Slope'};...
    'loom_spd_side_blk',...
        {'AnovaP';'MaxLoc';'SDI'};...
    'loom_spd_overh_blk',{};...
    'loom_spd_overh_wht',{};... 
    'GRATING_ORIENTATION',...
        {'DDI'};...     
    'DIRECTION_TUNING',...
        {'DirDI'};... 
    'OPTIC_FLOW_SPEED',... 
        {'AnovaP';'MaxLoc';'SDI'};...
    'OPTIC_FLOW_COHER',... 
        {'Slope';'PRaw'};...
};

%%% areas_to_show lists areas to be compared along the response parameters
%%% listed in vars_to_show
areas_to_show = {'por', 'p', 'li', 'lm'}; 

% tuning_summary_listfile contains protocol names in the first column and the
% corresponding summary filename of the file containing results for all analyzed 
% cells on this protocol in the second column
%%% recording_data_listfile contains info about recordings, including date
%%% (to use as subject ID) and recording location; first line is variable
%%% names with data below
tuning_summary_listfile = 'C:\Users\AM\Documents\Lab\recordings\tuning_summary_list.xlsx';
recordings_summary_filename = 'C:\Users\AM\Documents\Lab\Recordings Summary.xlsx'; % spreadsheet with info for each subject
cell_recordings_filename = 'C:\Users\AM\Documents\Lab\Cell Recordings'; % spreadsheet with info for each individual cell
spike_timing_file = 'C:\Users\AM\Documents\Lab\recordings\spike_timing_summary.mat'; % record of peristimulus spike timing info
runsByType_file = 'C:\Users\AM\Documents\Lab\recordings\runs_by_type.xlsx';
loom_type_labels = {'loom_spd_side_blk','loom_spd_overh_blk','loom_spd_overh_wht'}; % must match labels in cell_recordings file
signf_thresh = 0.05; % for labeling cells as significantly tuned by Anova p-value 

%% setup
vars_to_show = cell2table(vars_to_show(:,2)','VariableNames',vars_to_show(:,1)'); % make vars_to_show into table
for varind = 1:size(vars_to_show,2)
    if ischar(vars_to_show{1,varind}{:,:}) % if var is char rather than cell (because it was created with one cell element)
        vars_to_show{1,varind}{:,:} = {vars_to_show{1,varind}{:,:}}; % turn it from char into cell
    end
end

% get file lists if they do not exist yet
if ~exist('allfiles','var') || ~exist('files_with_dirs','var')
    fprintf('File lists not in workspace, running generate_filelist.m.\n')
    [allfiles files_with_dirs] = generate_filelist;
end 

% get recdate and cell name from each run - all IDs as numbers, not strings
% assumes files will have names of form 'm%Xc%Yr%Z', where %X = subject ID (date),
%   %Y = cell ID, %Z = run ID
cell_ids = inf(size(allfiles));
subject_ids = inf(size(allfiles));
run_ids = inf(size(allfiles));
for ind = 1:length(allfiles)
    thisfname = getfname(allfiles{ind});
    mloc = strfind(thisfname,'m');
    rloc = strfind(thisfname,'r');
    cloc = strfind(thisfname,'c');
    subject_ids(ind) = str2double(thisfname(mloc+1:cloc-1));
    cell_ids(ind) = str2double(thisfname(cloc+1:rloc-1));
    run_ids(ind) = str2double(thisfname(rloc+1:end));
end
rundata = table(subject_ids,cell_ids,run_ids,files_with_dirs,...
    'VariableNames',{'recdate','cell','run','path'});

%% put data from all analysis summaries into a table
[~, summaryFileNames] = xlsread(tuning_summary_listfile); 
summarydata = table(summaryFileNames(:,1),summaryFileNames(:,2),cell(size(summaryFileNames,1),1),...
    'VariableNames',{'protocol','summary_file','data'});
for ind = 1:height(summarydata)
    %%% import and format summary data for this protocol
    % use importdata because tdfread generates problems when there are extra
    % tabs at ends of lines
    importSummary = importdata(summarydata.summary_file{ind});
    varnames = importSummary.textdata(1,:); % variable names will be on the first line
    if size(varnames,2) == 1    % if all variable names got put into a single cell
        varnames = strsplit(varnames{:}); % split varnames into the separate variable names
    end
    varnames = regexprep(varnames,'^ ',''); varnames = regexprep(varnames,'^ ',''); % eliminate leading whitespace
    varnames = regexprep(varnames,' $',''); varnames = regexprep(varnames,' $',''); % eliminate trailing whitespace
    nonBlankVars = ~strcmp(varnames,'');
    varnames = varnames(nonBlankVars); % delete blank variables
    for varInd = 1:length(varnames)
        thisVarInstances = find(strcmp(varnames(varInd),varnames)); 
        if length(thisVarInstances) > 1 % if there are multiple vars with this name
            varnames{thisVarInstances(2)} = [varnames{varInd} '_2']; % rename to avoid duplicate table variables
        end
    end
    for varInd = 1:length(varnames)
        varnames(varInd) = strrep(varnames(varInd),'+','_plus_'); % replace invalid chars with variable-compatible chars
    end
    thisProtSummary = table(importSummary.textdata(2:end,1),'VariableNames',varnames(1)); % add file name strings from first column to table
    thisProtSummary = [thisProtSummary array2table(importSummary.data,'VariableNames',varnames(2:end))];
    summarydata.data(ind) = {thisProtSummary};
end

%% import recording data to label cell results
[~, ~, recs_summary] = xlsread(recordings_summary_filename); 
vartitles = recs_summary(1,:);
emptytitles = cell2mat(cellfun(@(x)~ischar(x),vartitles,'UniformOutput',false)); % check for blank variable names
vartitles(emptytitles) = {'blankvar'}; % replace blank names with 'blankvar'
vartitles = cellfun(@(x) strrep(x,' ','_'),vartitles,'UniformOutput',false); % make into table-compatible varnames
vartitles = cellfun(@(x) strrep(x,',','_'),vartitles,'UniformOutput',false); % make into table-compatible varnames
vartitles = cellfun(@(x) strrep(x,'.','_'),vartitles,'UniformOutput',false); % make into table-compatible varnames
vartitles = cellfun(@(x) strrep(x,'-','_'),vartitles,'UniformOutput',false); % make into table-compatible varnames
vartitles = cellfun(@(x) strrep(x,'(','_'),vartitles,'UniformOutput',false); % make into table-compatible varnames
vartitles = cellfun(@(x) strrep(x,')','_'),vartitles,'UniformOutput',false); % make into table-compatible varnames
recs_summary = cell2table(recs_summary(2:end,:),'VariableNames',vartitles);
emptyrows = cell2mat(cellfun(@(x)~ischar(x),recs_summary.date,'UniformOutput',false)); % check for blank rows
recs_summary(emptyrows,:) = []; % delete blank rows

% import data on individual cells, identify variable rows/columns
[~, ~, cell_recs] = xlsread(cell_recordings_filename);
cell_recs_num = cell_recs; % cell_recs with strings replaced by NaN for easier reading later
cell_recs_num(~cellfun(@isnumeric,cell_recs_num)) = {NaN}; % replace non-numeric cells with NaN
cell_recs_num = cell2mat(cell_recs_num); % convert to matrix
[subRowStart, subCol] = find(strcmp(cell_recs,'Subject_Name'));
[~, cellCol] = find(strcmp(cell_recs,'Cell_Num'));
[~, depthCol] = find(strcmp(cell_recs,'Depth_Measured_um')); % cell rf location relative to vert meridian
[~, areaCol] = find(strcmp(cell_recs,'Area')); % area from which this cell was recorded
[~, rfAzmCol] = find(strcmp(cell_recs,'Azim_Ctr')); % cell rf location relative to vert meridian
[~, rfElvCol] = find(strcmp(cell_recs,'Elev_Ctr')); % cell rf location relative to hrz meridian
[~, ofsAzmCol] = find(strcmp(cell_recs,'Azim_Ofs')); % screen azimuth offset from the coronal plane
[~, ofsElvCol] = find(strcmp(cell_recs,'Elev_Ofs')); % screen elevation offset from the coronal plane


% reformat dates to make compatible with subject names/dates in celldata
for ind = 1:height(recs_summary)
    recs_summary.date(ind) = strrep(recs_summary.date(ind),'/','');
    thisdate = recs_summary.date{ind}; 
    if length(thisdate) == 6; % if month and day are both 1 digit
        recs_summary.date{ind} = [thisdate(1),'0',thisdate(2:end)]; % add 0 before day
    elseif length(thisdate) == 7 && strcmp(thisdate(1),'1') % if month is 2 digits and day is 1 digit
        recs_summary.date{ind} = [thisdate(1:2),'0',thisdate(3:end)]; % add 0 before day
    end
    recs_summary.date{ind} = recs_summary.date{ind}([1:end-4,end-1:end]); % eliminate '20' on e.g. '2016'
    recs_summary.date{ind} = str2double(recs_summary.date{ind}); % make double
end

%% create a table row for each cell, including all summary data for this cell
celldata = table;
unq_subjs = unique(rundata.recdate,'stable');
for subind = 1:length(unq_subjs)
    thissubj=unq_subjs(subind);
    this_subj_rows = rundata(rundata.recdate==thissubj,:);
    unq_cells_thissubj = unique(this_subj_rows.cell);
    recdata_matchrow = find(cell2mat(recs_summary.date) == thissubj);  % find corresponding subject entry in recording_data
    for cellind = 1:length(unq_cells_thissubj)
        thiscell = unq_cells_thissubj(cellind);
        celldata = [celldata; table(thissubj,thiscell,{'_'},{'_'},NaN,NaN,NaN,NaN,NaN,'VariableNames',{'recdate','cell','subject','area','depth','rf_azm','rf_elv','ofs_azm','ofs_elv'})];
        % find matching data for this subject from recording_data
        if ~isempty(recdata_matchrow) % if corresponding entry was found in recording_data
            lastrow = height(celldata);
            subjID = recs_summary.name{recdata_matchrow};
            celldata.subject{lastrow} = subjID;
% % % % % % % %             celldata.area{lastrow} = recs_summary.probe_location{recdata_matchrow};     
        end
        % get data from the cell_recordings_file for this cell
        thisCellRow = find((cell_recs_num(:,cellCol) == thiscell) & (cell_recs_num(:,subCol) == subjID)); % row within cell_recs
        if ~isempty(thisCellRow) % if corresponding entry was found in cell recordings file, import info for this cell
            celldata.depth(lastrow) = cell_recs_num(thisCellRow,depthCol);
            celldata.area{lastrow} = cell_recs{thisCellRow,areaCol}; 
            celldata.rf_azm(lastrow) = cell_recs_num(thisCellRow,rfAzmCol); % RF center azimuth           
            celldata.rf_elv(lastrow) = cell_recs_num(thisCellRow,rfElvCol); % RF center elevation
            celldata.ofs_azm(lastrow) = cell_recs_num(thisCellRow,ofsAzmCol); % screen offset azimuth           
            celldata.ofs_elv(lastrow) = cell_recs_num(thisCellRow,ofsElvCol); % screen offset elevation
        end
    end
end

%% assign analysis results to corresponding cells
% protocol names must be found on second line of the cell_recordings file
ncells = height(celldata);
protNamesList = vars_to_show.Properties.VariableNames;
for indprottype = 1:length(protNamesList) % get list of protocol type each run is
    cols(indprottype) = find(strcmp(protNamesList(indprottype),cell_recs(2,:)));
end
prot_type_list = cell_recs(3:ncells+2,cols);
nonnumeric = find(cell2mat(cellfun(@(x)~isnumeric(x),prot_type_list,'UniformOutput',false))); 
prot_type_list(nonnumeric) = {NaN}; 
prot_type_list = cell2mat(prot_type_list);
prot_type_list = [celldata(:,'recdate'), celldata(:,'cell'), array2table(prot_type_list,'VariableNames',protNamesList)];

for protind = 1:height(summarydata)
    thisprot = summarydata.protocol{protind};
    thisprotdata = summarydata.data{protind};
    % make new variables in celldata
    if strcmp(thisprot,'LOOMING_SPEED') % add variables to celldata for all looming types
        celldata{:,loom_type_labels} = cell(height(celldata),length(loom_type_labels));        
    else
        celldata{:,thisprot} = cell(height(celldata),1); % add variable to celldata for this analysis protocol
    end
    % enter data from each run into celldata
    for protrowind = 1:height(thisprotdata)
        thisfname = thisprotdata.FILE{protrowind};
        mloc = strfind(thisfname,'m');
        cloc = strfind(thisfname,'c');
        rloc = strfind(thisfname,'r');
        dotloc = strfind(thisfname,'.htb');
        thissubj = str2double(thisfname(mloc+1:cloc-1));
        thiscell = str2double(thisfname(cloc+1:rloc-1));
        thisrun = str2double(thisfname(rloc+1:dotloc-1));
        matchrow = find(thissubj == prot_type_list.recdate &...
                        thiscell == prot_type_list.cell);
        if length(matchrow) > 1
            error('Multiple cell entries found in celldata for recdate %g, cell %g.',thissubj,thiscell)
        elseif length(matchrow) == 1
            colmatch = find(thisrun == prot_type_list{matchrow, protNamesList}); % find protocol for this run
            if length(colmatch) > 1 
                error('Multiple run entries found in celldata for recdate %g, cell %g, run %g.',thissubj,thiscell,thisrun)
            elseif length(colmatch) == 1
                thisprottitle = protNamesList{colmatch};
                celldata_row = find(celldata.recdate==thissubj & celldata.cell==thiscell);
                celldata{celldata_row,thisprottitle} = {thisprotdata(protrowind,:)};    
            end
        end
    end
end

%% tabulate variables of interest from vars_to_show for each cell

protlist = protNamesList; % protocol list
% % add variables to celldata
paramsToShowNames = {};
for protInd = 1:length(protlist)
	thisprot = protlist{protInd};    
    paramlist = vars_to_show{1,thisprot}{:}; % parameters of interest within this protocol
    for paramInd = 1:length(paramlist)
        thisparam = paramlist{paramInd};
        paramFullName = [thisprot,'_',thisparam];
        paramsToShowNames = [paramsToShowNames paramFullName]; % add param name to list of all parameters of interest
        celldata = [celldata table(NaN(ncells,1),'VariableNames',{[thisprot,'_',thisparam]})];   % add variable to celldata 
    end
end
% fill in data for individual cells
for row = 1:ncells
    for protInd = 1:length(protlist)
        thisprot = protlist{protInd}; 
        if ~isempty(celldata{row,thisprot}{:}) % if this protocol was run for this cell
            paramlist = vars_to_show{1,thisprot}{:}; % parameters of interest within this protocol
            thisSummary = celldata{row,thisprot}{:,:};
            for paramInd = 1:length(paramlist)
                thisparam = paramlist{paramInd};
                paramFullName = [thisprot,'_',thisparam];
                celldata{row,paramFullName} = thisSummary{1,thisparam}; % fill value into celldata
            end
        end
    end
end

% put summary table columns last within celldata
summaryVarNames = [summaryFileNames(:,1); protNamesList'];
summaryCols = cellfun(@(x)any(strcmp(x,summaryVarNames)),celldata.Properties.VariableNames); 
celldata = [celldata(:,~summaryCols) celldata(:,summaryCols)];

%% add spike timing info for each cell
% 1 indicates significantly different spike rate during response than
% during baseline, = nonsgn difference
load(spike_timing_file,'spikesTable');
rspvVarNames = cellfun(@(x)['rspv_',x],protNamesList,'UniformOutput',false);
responsiveTable = array2table(NaN(ncells,size(vars_to_show,2)),'RowNames',celldata.Properties.RowNames,'VariableNames',rspvVarNames);
celldata = [celldata responsiveTable]; 
for spikesTable_row = 1:height(spikesTable)
    thisfname = spikesTable.file_name{spikesTable_row};
    mloc = strfind(thisfname,'m');
    cloc = strfind(thisfname,'c');
    rloc = strfind(thisfname,'r');
    dotloc = strfind(thisfname,'.htb');
    thisdate = str2double(thisfname(mloc+1:cloc-1));
    thiscell = str2double(thisfname(cloc+1:rloc-1));
    thisrun = str2double(thisfname(rloc+1:dotloc-1));
    celldata_matchrow = find(thisdate == celldata.recdate(:) &...
                             thiscell == celldata.cell(:));
    if length(celldata_matchrow) > 1
        error('More than one match found in celldata for recdate %g, cell %g.',thissubj,thiscell)
    elseif length(celldata_matchrow) == 1 % if one matching row is found in celldata
        protTypeList_matchrow = find(thisdate == prot_type_list.recdate(:) &...
                                     thiscell == prot_type_list.cell(:));
        if length(protTypeList_matchrow) > 1 
            error('More than one match found in prot_type_list for recdate %g, cell %g.',thissubj,thiscell)
        elseif length(protTypeList_matchrow) == 1 % if one matching row is found in prot_type_list
            matchingRow = prot_type_list{protTypeList_matchrow,protNamesList};
            matchCol = find(thisrun == matchingRow);
            if length(matchCol) > 1
                error('More than one match found in prot_type_list for recdate %g, cell %g.',thissubj,thiscell)
            elseif length(matchCol) == 1
                thisProt = protNamesList{matchCol};
                celldata{celldata_matchrow,['rspv_',thisProt]} = spikesTable.sgn_resp(spikesTable_row);
            end
        end  
    end
end
                    

    
%% area comparison analysis
% means for each parameter for each area
area_means = array2table(NaN(length(areas_to_show),length(paramsToShowNames)),'VariableNames',paramsToShowNames,'RowNames',areas_to_show);
for indArea = 1:length(areas_to_show)
    thisarea = areas_to_show{indArea};
    for indparam = 1:length(paramsToShowNames)
        thisparam = paramsToShowNames{indparam};
        matchrow = strcmp(areas_to_show{indArea},celldata.area);
        thisParamThisArea = celldata{matchrow,thisparam};
        area_means{thisarea,thisparam} = nanmean(thisParamThisArea);
    end
end

% perform anova between area for each variable of interest
area_anova_p = array2table(NaN(size(paramsToShowNames)),'VariableNames',paramsToShowNames);
for indparam = 1:length(paramsToShowNames)  
    thisparam = paramsToShowNames(indparam);
    anovagroups = [];
    anovanames = {};
    for indArea = 1:length(areas_to_show)
        thisarea = areas_to_show{indArea};
        matchrow = strcmp(celldata.area,thisarea);
        anovagroups = [anovagroups; celldata{matchrow, thisparam}];
        ncells = length(find(matchrow));
        anovanames = [anovanames; cellstr(repmat((thisarea),ncells,1))];
    end
    area_anova_p{1,thisparam} = anova1(anovagroups,anovanames,'off');
end
        

%%
if clear_workspace
    clearvars -except celldata files_with_dirs area_means area_anova_p
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

