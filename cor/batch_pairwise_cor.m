%%% pairwise correlations for multiple planes/sessions
%
%  [cortable, filetable] = batch_pairwise_cor(filetable_in,pars)
%
%%%%%%% filetable_in = master excel file listing tuning and patchdata files to analyze, or table variable
%%%%%%%%% updated 2020/8/15

function [cortable, filetable] = batch_pairwise_cor(filetable_in,pars)

pars = vardefault('pars',struct);
pars.layers_to_analyze = field_default(pars,'layers_to_analyze', [  ]);
% pars.days_to_analyze =  field_default(pars, 'days_to_analyze', {'2018-02-02', '2018-03-28', '2018-04-06'}); % han lab
% pars.days_to_analyze =  field_default(pars, 'days_to_analyze', {'2018-12-30', '2019-01-04', '2019-01-07', '2019-01-18', '2019-01-19', '2019-01-20', '2019-01-21'}); % 4th floor scpoe
show_plots = 0; 
pars.show_waitbar = 0; % don't show wait bar within each session; doesn't affect waitbar in this function

if ~istable(filetable_in) % if filetable_in is a filename rather than an already-processed table variable
    filetable = excel2table(filetable_in); % load master file list
    filetable = filetable(filetable.analyze_plane == 1,:); % keep only planes that have be marked for analysis
    filetable.analyze_plane = [];  
elseif istable(filetable_in)
    filetable = filetable_in;
end
nfiles = height(filetable);
pars.days_to_analyze =  field_default(pars, filetable.day);

%%% get roi patch location data, concatenate in single table
cortable = table;
wbar = waitbar(0,'Gathering pairwise correlation data...');
for ifile = 1:nfiles
    if any(pars.layers_to_analyze == filetable.layer(ifile)) || isempty(pars.layers_to_analyze)
        if any(strcmp(pars.days_to_analyze, filetable.day{ifile})) || isempty(pars.days_to_analyze)
            tuning_file = [filetable.directory{ifile}, filesep, filetable.tuning_file{ifile}];
            tuningdat = load(tuning_file,'tuningdat'); tuningdat = tuningdat.tuningdat;
            %  if there's nostim_onset or nostim_offset are missing, use dark onset and dark offset
            if ~any(strcmp(filetable.Properties.VariableNames,'nostim_onset'))
                scopetiming = load([filetable.directory{ifile}, filesep, filetable.scopetiming_file_dark{ifile}]); % get dark period timing info
                nostim_onset = scopetiming.scopetiming.scope_start_timepoint;
            else
                nostim_onset = filetable.nostim_onset(ifile);
            end
            if ~any(strcmp(filetable.Properties.VariableNames,'nostim_offset'))
                scopetiming = load([filetable.directory{ifile}, filesep, filetable.scopetiming_file_dark{ifile}]); % get dark period timing info
                nostim_offset = scopetiming.scopetiming.scope_stop_timepoint;
            else
                nostim_offset = filetable.nostim_offset(ifile);
            end
            partial_cortable = pairwise_cor(tuningdat, nostim_onset, nostim_offset,  show_plots, pars);
            partial_cortable = [repmat(filetable(ifile,{'sub','day','plane','layer'}),height(partial_cortable),1), partial_cortable]; % add session data to roi table            
            cortable = [cortable; partial_cortable];
            try waitbar(ifile/nfiles,wbar); end
            clear tuning_file   partial_cortable partial_cortable
        end
    end
end
try close(wbar); end

%%% stats tests
