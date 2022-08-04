%TUNINGSTIM_PREPARESTIM
% Create table for recording stim parameters and stim presentation notes.
%%% last updated 1/28/16 on msi
function [stimrec] = tuningstim_preparestim(stimpars,commonvars)

%% Make stimulus record table
% uniqueTrials = number of unique parameter sets (values of the parameter
% to be tested)
% Set the variable to be tested for each trial.
switch stimpars.tuningParameter
    case 'orient'
        uniqueTrials = stimpars.n_orients;
        orient_table = table(stimpars.angle_vals','VariableNames',{'orient'});
    case 'sf'
        uniqueTrials = stimpars.n_sfs;
        sf_table = table(stimpars.sf_vals','VariableNames',{'sf'});
    case 'tf'
        uniqueTrials = stimpars.n_tfs;
        tf_table = table(stimpars.tf_vals','VariableNames',{'tf'});
    case 'size'
        uniqueTrials = stimpars.n_diams;
        diam_table = table(stimpars.diam_vals','VariableNames',{'diam'});
    case 'rf'
        uniqueTrials = stimpars.rows * stimpars.columns;
    case 'ssmod'
        uniqueTrials = stimpars.n_diams * stimpars.n_orients * (stimpars.n_tfs + stimpars.n_sfs); % may eliminate multiple orients for ssmod
end

singleVarNames = {'sf_pix','tf_degperflip','diam_pix'};% ,'missedFlips','precedingIsi','stimTimeStart','stimTimeEnd'};
singleVarNans = repmat({NaN(uniqueTrials,1)},1,length(singleVarNames));

% Set the constant variables for each trial.
if ~exist('orient_table','var')
    orient_table = table(stimpars.orient*ones(uniqueTrials,1),'VariableNames',{'orient'});
end
if ~exist('sf_table','var')
    sf_table = table(stimpars.sf*ones(uniqueTrials,1),'VariableNames',{'sf'});
end
if ~exist('tf_table','var')
    tf_table = table(stimpars.tf*ones(uniqueTrials,1),'VariableNames',{'tf'});
end
if ~exist('diam_table','var')
    diam_table = table(stimpars.diam*ones(uniqueTrials,1),'VariableNames',{'diam'});
end

nframes = round(stimpars.dur_stim / commonvars.ifi);
frameDuration_table = table(NaN(uniqueTrials,nframes),'VariableNames',{'frameDuration'}); % for checking for missed flips
stimrec = [orient_table sf_table tf_table diam_table...
    table(singleVarNans{:},'VariableNames',singleVarNames) frameDuration_table];

% Convert units (degs to pix, cycles/deg to cycles/pixel, secs to frames). We are assuming that 
% aspect ratio of screen resolution = aspect ratio of physical screen size (square pixels), 
% and therefore that we can use the same the same conversion for the horz and vert axes.  
%%% Note: deg2pixels is only valid for lines bisected by the screen center. 
%%%% degperflip uses deg along a sine wave oscillation, not deg in visual space
stimrec.sf_pix = 1./deg2pixels(1./stimrec.sf,commonvars); % convert degrees to pixels
stimrec.diam_pix = deg2pixels(stimrec.diam,commonvars); % convert degrees to pixels
stimrec.tf_degperflip = 360 * commonvars.ifi * stimrec.tf; %360deg/cyc * sec/flip *cyc/sec = deg/flip

stimrec = repmat(stimrec,stimpars.repetitions,1); %%% copy each unique param set stimpars.repetitions times
ntrials = height(stimrec);
stimrec = stimrec(randperm(ntrials),:); %%% shuffle trial order



    
    
    
