%WARPED_RUNSTIM - called from stim_main_warped
% Use WARPED_STIM_GENERATE beforehand to generate warped surround stimuli. 
%%% last updated 2019-1-18 on thermaltake

%% Setup
if ~exist('prerender_file','var')
    prerender_file = uigetfile('.mat','Select prerendered gratings file');
end
keepFrameDuration = 0; %%% save into the stim file the measured duration of every displayed frame
if save_stim_file
    stimpars.savefile = strcat('stimdata_',strrep(datestr(now),':',';'),'_tuning'); % name of stim log file
else
    checksave = input('Not saving stim file. Press Enter to continue.','s');
    stimpars.savefile = '';
end
load(prerender_file,'prerend')
setupvars.stimScreen = prerend.stimScreen;
commonvars = pulse_and_stim_setup(setupvars, prerender_file);
winPtr = commonvars.winPtr;
ifi = commonvars.ifi;
try stimpars = rmfield(stimpars,'stim_center'); end
stimpars.stim_center = flipud(fliplr(stimpars.stim_center_yx)); % psychtoolbox uses xy, so fliplr or flipud
% create and/or load stimulus gratings 
warped_preparestim; 

stimpars.sf_vals = sf_vals;
stimpars.tf_vals = tf_vals;
stimpars.orient_vals = orient_vals;
stimpars.diam_vals = diam_vals; 
stimpars.warping = 'warped'; % indicate that warped stimuli were used
stimpars.experiment_type = 'sf_tf_orient_diam'; 
stimpars.screenstats = commonvars.screenstats;
if commonvars.send_pulse
    stimpars.pulse = rmfield(pulsepars,{'stimon_signal','stimoff_signal'});
end
stimpars.computer = getenv('computername');
stimpars.prerender_file = prerender_file;
if commonvars.overwrite_check && (exist(stimpars.savefile,'file') || exist(strcat(stimpars.savefile,'.mat'),'file'))
    overwrite = input(sprintf('File ''%s'' already exists. Enter ''y'' to overwrite.\n', stimpars.savefile),'s');
    if ~strcmp(overwrite,'y')
        error('Will not overwrite file.')
    end
end
if save_stim_file
    save(stimpars.savefile,'stimpars','stimpar_sets');
end

if stimpars.grats_in_workspace
    fprintf('Loading prerendered warped gratings...')
    gratingLoadTic = tic;
    load(prerender_file,'outGratCells');
    fprintf([' took ' num2str(toc(gratingLoadTic)) 's.\n']);
end

stimpars.starttime = datestr(now);
warped_stimloop; %% present stimuli

stimpar_sets.precedingIsi(1) = NaN; % first trial has no preceding ISI

stimpar_sets.grating_table_rows = []; % clear unnecessary variables from stimpar_sets
stimpar_sets.outer_apertex = [];
stimpar_sets.outer_apt_rect = [];
if stimpars.SHOW_CENTER
    stimpar_sets.inner_gratingtex = [];
end
if save_stim_file
    stimpar_sets.outerGratSrcRect = [];
    if ~keepFrameDuration 
        stimpar_sets.frameDuration = [];
    end
    save(stimpars.savefile,'stimpars','stimpar_sets')
    disp(sprintf(['Stimulus presentation complete.',...
        '\nStimulus data saved in ' stimpars.savefile '.mat.'...
        '\nEnd recording on acquisition computer.']));
end

% Cleanup after stim presentation
Screen('CloseAll')
Screen('Preference','SkipSyncTests',commonvars.origSyncLevel); % reset to original sync testing setting
Screen('Preference', 'VisualDebugLevel', commonvars.origDebugLevel); % reset to iriginal PTB debug level
Screen('Preference','Verbosity',commonvars.origVerbosity); % reset to original PTB warning mode
commandwindow;
clear global