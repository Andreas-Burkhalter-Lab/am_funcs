%%%%%%%%%% load data from files for tuning curve processing
% part of get_tuning_curves
% updated 5/29/18 on thermaltake

%% load files, set params
if ~exist('dff_file','var') || isempty(dff_file)
    dff_file = uigetfile('*.mat','Select dff .mat file.','MultiSelect','off');
end

if ~exist('abf_file','var') || isempty(abf_file)
    abf_file = uigetfile('*.abf','Select .abf file.');
end

if ~exist('stimdata_file','var') || isempty(stimdata_file)
    stimdata_file = uigetfile('*stimdata*.mat','Select stimdata .mat file.');
end

if ~exist('scopetiming_file','var') || isempty(scopetiming_file)
    scopetiming_file = uigetfile('*_scopetiming.mat','Select scopetiming .mat file.');
end

if ~exist('regops_file') || isempty(regops_file)
    regops_file = uigetfile('*regops*.mat','Select regops .mat file.');
end

try
    [abf_timepoints abf_sampling_interval_us abfinfo] = abfload(abf_file); % abf_sampling_interval_us = sampling interval microseconds
catch abf_mexception
    fprintf('     ABF file not readable as .abf; loading as .mat file instead.\n')
    load(abf_file,'-mat')
end
load(stimdata_file);
load(scopetiming_file);
load(regops_file); 
dff_data = load(dff_file,'dff','meanImage','masks','fName');
if exist([getfname(dff_data.fName) '.mat'],'file')
    planeinfo = load([getfname(dff_data.fName) '.mat'],'ops');
    iplane = planeinfo.ops.iplane; % plane number for aligning with TTL pulses; if iplane==1, first frame might be wrong
else 
    iplane = input('Enter plane number: ');
end

if isfield(ops1{1},'nplanes')
    nplanes_in_session_total = ops1{1}.nplanes;
else 
    nplanes_in_session_total = ops1{1}.numPlanes;
end
stimdur_us = pars.stimdur * 1e6; % convert from s to microseconds
stimdur_scans = round(stimdur_us / abf_sampling_interval_us);
ntrials = height(par_sets);
nimages = size(dff_data.dff,1);
ncells = size(dff_data.dff,2);
res.stimpars = pars;
res.iplane = iplane;
res.nplanes_in_session_total = nplanes_in_session_total;

% get stim events from abf file
stimchan = find(strcmp(abfinfo.recChNames,'Psycho_In'));
locm_forw_chan = find(strcmp(abfinfo.recChNames,'Ball_forw'));
locm_yaw_chan = find(strcmp(abfinfo.recChNames,'Ball_yaw'));
stim_events = find_stimchan_events(abf_timepoints(:,stimchan)',pars,par_sets,abf_sampling_interval_us);
stim_events.Properties.VariableNames{find(strcmp('onset',stim_events.Properties.VariableNames))} = 'stim_onset';
stim_events.Properties.VariableNames{find(strcmp('duration',stim_events.Properties.VariableNames))} = 'stim_duration';
stim_events.chan_val = [];
par_sets = [par_sets stim_events]; clear stim_events
