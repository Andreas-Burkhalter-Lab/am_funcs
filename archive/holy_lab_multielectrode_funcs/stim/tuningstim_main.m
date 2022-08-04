%TUNINGSTIM_MAIN % Present stimuli to test tuning preferences.  
% Present circular gratings of with one parameter varying; parameter can be
% orientation, spatial frequency, temporal frequency, or stimulus size.
% A PTB window in which stimuli are to be presented should be open when
%       this script is called.  
% .... Gao et al. had 2-degree-diameter patches or gratings, 2s stim
%      presentation
%%% Last updated 1/28/16
function [stimrec, stimpars_out, stimpref_out] =...
    tuningstim_main(stimpars_in,stimpref_in,pulsepars,defaultpars,commonvars)

%% Setup
stimpars_out = stimpars_in; stimpref_out = stimpref_in;
global daq_sess
stimpars_out.warping = 'nonwarped'; % indicate nonwarped stimuli used for orientation tuning testing
singleVStimOn = 1; % only one vstimon value will be used

%%% assign values for the parameter being tested to present
switch stimpars_out.tuningParameter
    case 'orient'
        angle_increment = 360/stimpars_out.n_orients; % degrees of spacing between angles to present
        stimpars_out.angle_vals = linspace(stimpars_out.zeroorient,360-angle_increment,stimpars_out.n_orients);
    case 'sf'
        log_sf_minmax = log(stimpars_out.sf_minmax);
        log_sf_range = linspace(log_sf_minmax(1),log_sf_minmax(2),stimpars_out.n_sfs);
        stimpars_out.sf_vals = exp(log_sf_range);
    case 'tf'
        log_tf_minmax = log(stimpars_out.tf_minmax);
        log_tf_range = linspace(log_tf_minmax(1),log_tf_minmax(2),stimpars_out.n_tfs);
        stimpars_out.tf_vals = exp(log_tf_range);
    case 'diam'
        log_diam_minmax = log(stimpars_out.diam_minmax);
        log_diam_range = linspace(log_diam_minmax(1),log_diam_minmax(2),stimpars_out.n_diams);
        stimpars_out.diam_vals = exp(log_diam_range);
end
        
% Set stim pars to defaults and stim preferences already measured.
stimpars_out = copyAllFields(stimpars_out,defaultpars); % use default parameters for all variables except the one being tested
stimpars_out = copyAllFields(stimpars_out,stimpref_out); % use tuning preferences which have already been tested

% run stim checks and make stim record table
pulsepars = check_stim_get_pulse(stimpars_out,pulsepars,daq_sess,singleVStimOn,commonvars);
stimrec = prepare_tuning_stim(stimpars_out,commonvars); 

% Save log file
stimpars_out.logfile = [stimpars_out.tuningParameter 'tuning_stimlog_' date];
if commonvars.overwrite_check && (exist(stimpars_out.logfile,'file') || exist(strcat(stimpars_out.logfile,'.mat'),'file'))
    overwrite = input(sprintf('File ''%s'' already exists. Enter ''y'' to overwrite.\n', stimpars_out.logfile),'s');
    if ~strcmp(overwrite,'y')
        error('Will not overwrite file.')
    end
end
stimpars = stimpars_out; % rename for logfile
stimrec = table2dataset(stimrec); % matlab 2012 can handle datasets, not tables
save(stimpars_out.logfile,'stimrec','stimpars','pulsepars','commonvars');
stimrec = dataset2table(stimrec); 

%% Present stim
stimpars_out.ntrials = height(stimrec);
presentation_dur_minutes = (stimpars_out.dur_stim+stimpars_out.isi) * stimpars_out.ntrials * stimpars_out.repetitions / 60;
input('Begin Merec2 recording on diesel, then press Enter on stimulus computer.');
fprintf(['\nNow presenting %s tuning stimuli...\n',...
    'Approximate duration = %g minutes.\n'],...
    stimpars_out.tuningParameter,presentation_dur_minutes);

tuningstim_stimloop(stimrec,stimpars_out,pulsepars,commonvars); % stim presentation loop

save(stimparsout.logfile,'stimPresentationStartTime','-append')


%% Have user input tuning preference determined on recording comp. 
switch stimpars_out.tuningParameter
    case 'orient'
        stimpref_out.orient = input('Enter preferred orientation (deg from upright) calculated on recording computer:\n'); 
    case 'sf'
        stimpref_out.sf = input('Enter preferred spatial frequency (cyc/deg) calculated on recording computer:\n'); 
    case 'tf'
        stimpref_out.tf = input('Enter preferred temporal frequency (hz) calculated on recording computer:\n'); 
    case 'diam'
        stimpref_out.diam = input('Enter preferred stimulus diameter (deg) calculated on recording computer:\n'); 
end






