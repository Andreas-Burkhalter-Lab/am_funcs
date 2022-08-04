%%% this version separates out and plots variables of interest differently
%%% from newer version

% create a table with a row containing all data relevant to a given
% cell; run this script after generate_filelist and getTempoTrialData
%%% 8/2/16 AM edited on msi

%% parameters
%%% vars_to_show lists parameters of interest to be presented in summary
%%% tables
% first column in vars_to_show is the protocol, second is parameter names
% within protocol summary files
vars_to_show = {
    'RF_MAPPING',{};...    %%% rf mapping with drifting dots, not gratings
    'GRATING_RF_MAP',{};... 
    'SIZE_TUNING',{};... %%% dot or dot field size, not grating size
    'GRATING_SIZE',{};... 
    'GRATING_SPATIAL_FREQ',...
        {'AnovaP';'MaxLoc';'SFDI'};...
    'GRATING_TEMPORAL_FREQ',...
        {'AnovaP';'MaxLoc',;'TFDI'};...
    'SPEED_TUNING',... % drifting dot field speed
        {'AnovaP';'MaxLoc',;'SDI'};...
    'MOTION_COHERENCE',... % drifting dot field motion coherence
        {'PRaw';'Slope'};...
    'LOOMING_SPEED',{};... 
    'GRATING_ORIENTATION',{};... 
    'DIRECTION_TUNING' ,{};... 
    'OPTIC_FLOW_SPEED',{};... 
    'OPTIC_FLOW_COHER',{};... 
};

%%% area_to_show lists areas to be compared along the response parameters
%%% listed in vars_to_show
areas_to_show = {'por', 'p', 'li'}; 

% summary_file_listfile contains protocol names in the first column and the
% corresponding summary filename of the file containing results for all analyzed 
% cells on this protocol in the second column
%%% recording_data_listfile contains info about recordings, including date
%%% (to use as subject ID) and recording location; first line is variable
%%% names with data below
summary_file_listfile = 'C:\Users\AM\Documents\Lab\recordings\summary_file_list';
recordings_summary_filename = 'C:\Users\AM\Documents\Lab\Recordings Summary.xlsx'; % spreadsheet with info for each subject
cell_recordings_filename = 'C:\Users\AM\Documents\Lab\Cell Recordings'; % spreadsheet with info for each individual cell
runsByType_file = 'C:\Users\AM\Documents\Lab\recordings\runs_by_type.xlsx';
loom_type_labels = {'loom_spd_side_blk','loom_spd_overh_blk','loom_spd_overh_wht'}; % must match labels in runsByType_file
signf_thresh = 0.05; % for labeling cells as significantly tuned by Anova p-value 

%% setup
vars_to_show = cell2table(vars_to_show(:,2)','VariableNames',vars_to_show(:,1)'); % make vars_to_show into table

% get file lists if they do not exist yet
if ~exist('allfiles','var') || ~exist('files_with_dirs','var')
    fprintf('File lists not in workspace, running generate_filelist.m.\n')
    [allfiles files_with_dirs] = generate_filelist;
end 

% get recdate and cell name from each run - all IDs as numbers, not strings
% assumes files will have names of form 'm%Xc%Yr%Z', where %X = ` ID (date),
%   Y = cell ID, %Z = run ID
cell_ids = inf(size(allfiles));
subject_ids = inf(size(allfiles));
run_ids = inf(size(allfiles));
for ind = 1:length(allfiles)
    thisfname = getfname(allfiles{ind});
    mloc = strfind(thisfname,'m');
    cloc = strfind(thisfname,'c');
    rloc = strfind(thisfname,'r');
    subject_ids(ind) = str2double(thisfname(mloc+1:cloc-1));
    cell_ids(ind) = str2double(thisfname(cloc+1:rloc-1));
    run_ids(ind) = str2double(thisfname(rloc+1:end));
end
rundata = table(subject_ids,cell_ids,run_ids,files_with_dirs,...
    'VariableNames',{'recdate','cell','run','path'});

%% put data from all analysis summaries into a table
[~, summaryFileNames] = xlsread(summary_file_listfile); 
summarydata = table(summaryFileNames(:,1),summaryFileNames(:,2),cell(size(summaryFileNames,1),1),...
    'VariableNames',{'protocol','summary_file','data'});
for ind = 1:height(summarydata)
    %%% immport and format summary data for this protocol
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
[~, rfAzmCol] = find(strcmp(cell_recs,'Azim_Ctr')); % cell rf location relative to vert meridian
[~, rfElvCol] = find(strcmp(cell_recs,'Elev_Ctr')); % cell rf location relative to hrz meridian
[~, OfsAzmCol] = find(strcmp(cell_recs,'Azim_Ofs')); % screen azimuth offset from the coronal plane
[~, OfsElvCol] = find(strcmp(cell_recs,'Elev_Ofs')); % screen elevation offset from the coronal plane

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
unq_subjs = unique(rundata.recdate);
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
            celldata.area{lastrow} = recs_summary.probe_location{recdata_matchrow};       
        end
        % get data from the cell_recordings_file for this cell
        thisCellRow = find((cell_recs_num(:,cellCol) == thiscell) & (cell_recs_num(:,subCol) == subjID)); % row within cell_recs
        if ~isempty(thisCellRow) % if corresponding entry was found in cell recordings file, import info for this cell
            celldata.depth(lastrow) = cell_recs_num(thisCellRow,depthCol);
            celldata.rf_azm(lastrow) = cell_recs_num(thisCellRow,rfAzmCol); % RF center azimuth           
            celldata.rf_elv(lastrow) = cell_recs_num(thisCellRow,rfElvCol); % RF center elevation
            celldata.ofs_azm(lastrow) = cell_recs_num(thisCellRow,ofsAzmCol); % screen offset azimuth           
            celldata.ofs_elv(lastrow) = cell_recs_num(thisCellRow,ofsElvCol); % screen offset elevation
        end
    end
end

% sort looming files by type
[~, ~, runs_by_type] = xlsread(runsByType_file);
vartitles = runs_by_type(1,:);
runs_by_type = cell2table(runs_by_type(2:end,:),'VariableNames',vartitles);
loom_type_list = runs_by_type(:,loom_type_labels); 

% assign analysis results to corresponding cells
for protind = 1:height(summarydata)
    thisprot = summarydata.protocol{protind};
    thisprotdata = summarydata.data{protind};
    % make new variables in celldata
    if strcmp(thisprot,'LOOMING_SPEED') % add variables to celldata for peak, significance, and hwhm for all looming types
        celldata{:,loom_type_labels} = cell(height(celldata),length(loom_type_labels));
        loom_sgn_labels = cellfun(@(x)[x,'_sgnf'],loom_type_labels,'UniformOutput',false);
        celldata{:,loom_sgn_labels} = NaN(height(celldata),length(loom_sgn_labels));
        loom_peak_labels = cellfun(@(x)[x,'_peak'],loom_type_labels,'UniformOutput',false);
        celldata{:,loom_peak_labels} = NaN(height(celldata),length(loom_peak_labels)); 
        loom_hwhm_labels = cellfun(@(x)[x,'_hwhm'],loom_type_labels,'UniformOutput',false);
        celldata{:,loom_hwhm_labels} = NaN(height(celldata),length(loom_hwhm_labels));        
        hwhm_varind = find(~cellfun(@isempty,regexpi(thisprotdata.Properties.VariableNames,'hwhm')),1); % get hwhm variable index; use first matching variable
        hwhm_varname = thisprotdata.Properties.VariableNames{hwhm_varind}; 
    else
        celldata{:,thisprot} = cell(height(celldata),1); % add variable to celldata for this analysis protocol
        if any(find(strcmpi(thisprotdata.Properties.VariableNames,'AnovaP'))) % if AnovaP was computed for this protocol
            celldata{:,[thisprot '_sgnf']} = NaN(height(celldata),1); % add variable to celldata for this analysis protocol
        end
        peakvarind = find(strcmpi(thisprotdata.Properties.VariableNames,'peak'));
        if any(peakvarind) % if peak was computed for this protocol
            celldata{:,[thisprot '_peak']} = NaN(height(celldata),1); % add variable to celldata for this analysis protocol
            peakvarname = thisprotdata.Properties.VariableNames{peakvarind}; 
        end
        hwhm_varind = find(~cellfun(@isempty,regexpi(thisprotdata.Properties.VariableNames,'hwhm')),1); % get hwhm variable index; use first matching variable
        if any(hwhm_varind) % if hwhm was computed for this protocol
            celldata{:,[thisprot '_hwhm']} = NaN(height(celldata),1); % add variable to celldata for this analysis protocol
            hwhm_varname = thisprotdata.Properties.VariableNames{hwhm_varind}; 
        end
    end
    % enter data from each run into celldata
    for protrowind = 1:height(thisprotdata)
        thisfname = thisprotdata.FILE{protrowind};
        mloc = strfind(thisfname,'m');
        cloc = strfind(thisfname,'c');
        rloc = strfind(thisfname,'r');
        thissubj = str2double(thisfname(mloc+1:cloc-1));
        thiscell = str2double(thisfname(cloc+1:rloc-1));
        % if looming protocol, label it with the type of looming trial
        if strcmp(thisprot,'LOOMING_SPEED')
            matchtype = any(strcmp(thisfname,loom_type_list{:,:}));
            if any(matchtype)
                thisprottitle = loom_type_list.Properties.VariableNames{matchtype};
            end
        else
            thisprottitle = thisprot;
        end
        celldata_row = find(celldata.recdate==thissubj & celldata.cell==thiscell);
        celldata{celldata_row,thisprottitle} = {thisprotdata(protrowind,:)};
        if any(find(strcmp(thisprotdata.Properties.VariableNames,'AnovaP'))) % if AnovaP was computed for this protocol
            celldata{celldata_row,[thisprottitle '_sgnf']} = thisprotdata.AnovaP(protrowind) < signf_thresh;
        end    
        if any(peakvarind) % if peak was computed for this protocol
            celldata{celldata_row,[thisprottitle '_peak']} = thisprotdata{protrowind,peakvarname}; % list peaks for both tuned and untuned cells
        end    
        if any(hwhm_varind) % if half-width at half maximum was computed for this protocol
            celldata{celldata_row,[thisprottitle '_hwhm']} = thisprotdata{protrowind,hwhm_varname};
        end
    end
end

% make table showing proportion of significantly tuned cells for each stim parameter
sign_vars = ~cellfun(@isempty,(strfind(celldata.Properties.VariableNames,'_sgnf')));
sgn_table = celldata(:,sign_vars);
validcells = sum(~isnan(sgn_table{:,:}))'; % number of cells which were analyzed for a stim parameter 
n_sign_cells = nansum(sgn_table{:,:})'; % number of cells signficantly tuned for a stim param
fraction_sign = (n_sign_cells./validcells); % fraction of cells analyzed for a stim param that were signf tuned for it
sigf_summary = table(validcells,n_sign_cells,fraction_sign,'RowNames',sgn_table.Properties.VariableNames');

% reorganize table
peak_vars = ~cellfun(@isempty,(strfind(celldata.Properties.VariableNames,'_peak')));
other_vars = ~(sign_vars | peak_vars);
celldata = [celldata(:,peak_vars), celldata(:,sign_vars), celldata(:,other_vars)]; % significance and peak results before full data
leadingvars = strcmp(celldata.Properties.VariableNames,'subject') | strcmp(celldata.Properties.VariableNames,'recdate') |...
    strcmp(celldata.Properties.VariableNames,'cell') | strcmp(celldata.Properties.VariableNames,'area');
celldata = [celldata(:,{'subject' 'recdate' 'cell' 'area'}) celldata(:,~leadingvars)]; % most important variables first

% sort by area
indspor = find(strcmp(celldata{:,'area'},'por'));
indsli = find(strcmp(celldata{:,'area'},'li'));
indsother = [1:height(celldata)]'; indsother([indspor; indsli]) = [];
celldata = celldata([indspor; indsli; indsother],:); 


    
%% plotting
% subvariable, choose peak or hwhm
do_plotting = 0; 

if do_plotting
    close all
    tuned_cells_only = 1; % only plot cells which are significantly tuned
    % vartoplot = 'OPTIC_FLOW_SPEED';
    % vartoplot = 'loom_spd_side_blk';
    vartoplot = 'GRATING_SPATIAL_FREQ';
    % vartoplot = 'GRATING_TEMPORAL_FREQ';
    subvariable = 'peak';
    % subvariable = 'hwhm';
    tuning_var = [vartoplot '_sgnf'];
    tuned_cells = celldata{:,tuning_var}; 
    tuned_cells(isnan(tuned_cells)) = false;
    indspor = strcmp(celldata{:,'area'},'por');
    indsli = strcmp(celldata{:,'area'},'li');
    indsp = strcmp(celldata{:,'area'},'p');
    indsother = [1:height(celldata)]'; indsother([indspor | indsli | indsp]) = [];
    var_plus_subvar = [vartoplot '_' subvariable];
    if tuned_cells_only
        datpor = celldata{indspor&tuned_cells,var_plus_subvar};
        datli = celldata{indsli&tuned_cells,var_plus_subvar};
        datp = celldata{indsp&tuned_cells,var_plus_subvar};
    else 
        datpor = celldata{indspor,var_plus_subvar};
        datli = celldata{indsli,var_plus_subvar};
        datp = celldata{indsp,var_plus_subvar};
    end
    scatter([ones(size(datpor)); 2*ones(size(datli)); 3*ones(size(datp))], [datpor;datli;datp]);
    xlim([0.5 3.5]);
    title(var_plus_subvar);
    set(gca,'XTickLabel',{[] 'POR' [] 'LI' [] 'P' []}); 
    [h_por_vs_li pvalue_por_vs_li] = ttest2(datpor,datli);
    [h_por_vs_p pvalue_por_vs_p] = ttest2(datpor,datp);
end
    
    
    
    
    
    
    
    
    
    

