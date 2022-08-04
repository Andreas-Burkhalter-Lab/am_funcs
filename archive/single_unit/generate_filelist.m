% get .htb filenames from 'Cell Recordings' of the specified protocols into a list for getTempoTrialData
%%% last edited 8/13/16 on msi

function [files_with_dirs allfiles files_with_protocols] = generate_filelist

% topdir should contain folders containing .htb files 
topdir = 'C:\Users\AM\Documents\Lab\recordings'; 
cell_recordings_filename = 'C:\Users\AM\Documents\Lab\Cell Recordings'; % spreadsheet with info for each individual cell

% only the runs listed under the following protocol headings in 'Cell Recordings' will
% be analyzed; these names need to occur on the second line of the file
%%% do not include recede_spd_overh_blk because these cannot yet be
%%% read by looming speed tuning analysis script
protocols_to_analyze = {'GRATING_RF_MAP','SIZE_TUNING','GRATING_SPATIAL_FREQ','GRATING_TEMPORAL_FREQ','SPEED_TUNING','MOTION_COHERENCE','loom_spd_side_blk',...
    'GRATING_ORIENTATION','DIRECTION_TUNING','OPTIC_FLOW_SPEED','OPTIC_FLOW_COHER','loom_spd_overh_blk','loom_spd_overh_wht'};
[~, ~, cell_recs] = xlsread(cell_recordings_filename);
cell_recs = cell_recs(2:end,:); % eliminate first line; protocol titles should be on second line
extraVars = {'Date','Subject_Name','Cell_Num'}; % vars for finding file names/dirs
colHeadingsToKeep = [extraVars protocols_to_analyze];
keepCellRecsCols = [];
for cellRecsColumn = 1:size(cell_recs,2)
    if any(strcmp(cell_recs{1,cellRecsColumn},colHeadingsToKeep)) % if the column heading is a protocol name
        keepCellRecsCols = [keepCellRecsCols cellRecsColumn]; % list the column to keep
    end
end
cell_recs = cell_recs(:,keepCellRecsCols);
cell_recs = cell2table(cell_recs(2:end,:),'VariableNames',cell_recs(1,:)); % convert to table
datesAreStrings = cell2mat(cellfun(@isstr,cell_recs.Date,'UniformOutput',false)); % find valid rows of cell_recs
cell_recs = cell_recs(datesAreStrings,:); % delete extra rows
cell_recs.Date = cellstr(datestr(datenum(cell_recs.Date,'mm/dd/yyyy'),'mmddyy',2000)); % convert date format to that of directories
for protind = 1:length(protocols_to_analyze) % within run numbers, set non-numerics to NaN and convert all cells to matrices
    thisprot = protocols_to_analyze{protind};
    if iscell(cell_recs{:,thisprot})
        non_numerics = ~cellfun(@isnumeric,cell_recs{:,thisprot}); 
        cell_recs{non_numerics,thisprot} = {NaN}; % replace non-numeric cells with NaN
        cell_recs{non_numerics,thisprot} = {NaN};
        newvar = cell_recs{:,thisprot};
        cell_recs(:,thisprot) = []; % delete var
        cell_recs = [cell_recs table(cell2mat(newvar),'VariableNames',{thisprot})]; % re-add variable as mat with NaNs
        cell_recs = [cell_recs(:,1:length(extraVars)+protind-1) cell_recs(:,thisprot) cell_recs(:,length(extraVars)+protind:end-1)]; %reorer
    end
end
protocolRows = cell_recs{:,length(extraVars)+1:end}; 

% generate expected directory/filenames and check that they exist within
% topdir
nFiles = length(find(~isnan(protocolRows))); % number of files to analyze
allfiles = cell(nFiles,1);
files_with_dirs = cell(nFiles,1);
files_with_protocols = table(allfiles,allfiles,'VariableNames',{'filename','protocol'});

fileInd = 0;
for rowInd = 1:height(cell_recs)
    thisSub = cell_recs.Subject_Name(rowInd);
    for protInd = 1:length(protocols_to_analyze)
        thisprot = protocols_to_analyze{protInd};
        if ~isnan(cell_recs{rowInd,thisprot}) % if a valid run for this subject/protocol is specified
            fileInd = fileInd+1;
            thisdate = regexprep(cell_recs.Date{rowInd},'^0',''); % eliminate leading zero
            filename = sprintf('m%sc%gr%g.htb', thisdate, cell_recs.Cell_Num(rowInd), cell_recs{rowInd,thisprot}); 
            allfiles{fileInd} = filename;
            files_with_protocols.filename(fileInd) = {filename};
            files_with_protocols.protocol(fileInd) = {thisprot};
            files_with_dirs{fileInd} = sprintf('%s\\%g - %s\\%s',topdir,thisSub,cell_recs.Date{rowInd},filename);
            if ~exist(files_with_dirs{fileInd},'file')
                error(['Did not find recording ' files_with_dirs{fileInd} '.'])
            end
        end
    end
end























