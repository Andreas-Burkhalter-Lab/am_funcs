%----------------------------------------------------------------------------------------------------------
%-- AnalysisSwitchyard.m: This function takes in the data, protocol, etc. and chains off to the 
%-- 	relevant protocol-specific analysis code.  GCD, 1/3/2000
%--     Last Revised BJP, 3/1/02
%----------------------------------------------------------------------------------------------------------

%%% 5/19/16 AM added conditional in order to call different size tuning
%%% analyses depending on whether the protocol varied dot size or
%%% grating size
%%%% 8/29/16 AM added optional 'skip_trials' argin for GratingSF_Analyses/'Plot Tuning Curve' - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped
%%% 9/27/16 AM routed looming stimuli to ComputeSpikeRates_looming.m rather
%%%     than ComputeSpikeRates.m to change spike integration time based on
%%%     looming expansion rate; toggle this feature with loomingSpeedDependentSpikeWindow  

function Analysis_Switchyard(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffsetTime, StopOffsetTime, PATH, FILE, batch_flag, UseSyncPulses, skip_trials)

% AM added conditional to make skip_trials blank if unspecified)
if ~exist('skip_trials','var')
    skip_trials = [];
end

loomingSpeedDependentSpikeWindow = 1; % % AM 9/27/16; if on, spike counting window for looming stim will depend on looming speed
TEMPO_Defs;
ProtocolDefs;
Path_Defs;	%this sets up paths to the analysis functions, both Common and ProtocolSpecific

if(batch_flag == 1)
    ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
    t=sprintf('Doing analysis (%s) for Protocol: %s', Analysis{1}, protocol_names{Protocol+1});
    set(ListHandle, 'String', t);
end

%following function checks offset times and outputs starting and ending analysis (in spike bins)
num_trials = size(data.event_data, 3);
[StartOffsetBin StopOffsetBin StartEventBin StopEventBin] = CheckTimeOffset(data, num_trials, StartCode, StopCode, StartOffsetTime, StopOffsetTime, UseSyncPulses);
data.UseSyncPulses=UseSyncPulses;

%%% 9/27/16 AM routed looming stimuli to ComputeSpikeRates_looming rather than 
%%%   ComputeSpikeRates to change spike integration time based on looming expansion rate    
if (~isempty(data.spike_data))  
    %compute the firing rate over all trials during the period between StartCode and StopCode
    if loomingSpeedDependentSpikeWindow && any(Protocol == [LOOMING_SPEED; LOOMING_GAMP; LOOMING_GASSIAN_AMP]) % if it's a looming protocol and loomingSpeedDependentSpikeWindow==1
        data.spike_rates = ComputeSpikeRates_looming(data, num_trials, StartCode, StopCode, StartOffsetBin, StopOffsetBin); % use spike window based on looming speed        
    elseif ~loomingSpeedDependentSpikeWindow || ~any(Protocol == [LOOMING_SPEED; LOOMING_GAMP; LOOMING_GASSIAN_AMP]); % if it's not looming or we are using fixed spike windows for looming
        data.spike_rates = ComputeSpikeRates(data, num_trials, StartCode, StopCode, StartOffsetBin, StopOffsetBin); % use same spike window for all stim params
    end
end

if (~isempty(data.eye_data))
    %add the compute eye pos function here
    data.eye_positions = ComputeMeanEyePos(data, num_trials, StartCode, StopCode, StartOffsetBin, StopOffsetBin);
    %this next function checks to see if there is an eye calibration parameter file.  IF it exists,
    %the function will load in the parameters and compute the calibrated eye positions, and then store 
    %them back in the 'data' structure for use elsewhere.  Added 12/6/00 by GCD
    data.eye_calib_done = 0;  % set this flag to zero; it can be used elsewhere to check if there are calibrated signals
%    [caldata, doneflag] = LoadEyeCalibration(data, PATH, FILE);
%    [caldata, doneflag] = LoadEyeCalibration_NonLin(data, PATH, FILE);
%    data.eye_positions_calibrated = caldata;
%    data.eye_calib_done = doneflag;
end

switch(Protocol)		%call a .m file that contains protocol-specific analysis routines
case DIRECTION_TUNING
    DirectionTuning_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, skip_trials);
case SPEED_TUNING
    SpeedTuning_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, skip_trials);
case SIZE_TUNING
   
    %%% 5/19/16 AM added conditional in order to call different size tuning
    %%% analyses depending on whether the protocol varied dot size or
    %%% grating size; previously, SizeTuning_Analysis was called without
    %%% checking which stimulus type was used
    nDotSizes = length(unique(data.dots_params(DOTS_AP_XSIZ,:,PATCH1)));
    if nDotSizes > 1 % if dot size was varied
        DotSizeTuning_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE, skip_trials);
    else % if grating size was varied
        SizeTuning_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE, skip_trials); % previously only this line was included
    end   
    
case RF_MAPPING
    RFMapping_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, PATH, FILE, skip_trials);
case GRATING_ORIENTATION
    GratingOrientation_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, skip_trials);
case GRATING_SPATIAL_FREQ
    GratingSF_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, skip_trials);
case GRATING_TEMPORAL_FREQ
    GratingTF_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, skip_trials);
case GRATING_CONTRAST
    GratingCT_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
case GRATING_SIZE
    GratingSize_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, skip_trials);
case MOTION_COHERENCE
    MotionCoherence_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, skip_trials);
case GRATING_RF_MAP
    GratingRFMap_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, skip_trials);

%-------------------------------------------
% DoHyun Kim summer work at Dr. Burkhalter's Lab, WashU Medical School -
% Summer 2012 %fix later
%-------------------------------------------
case LOOMING_GAMP 
    LoomingContrast_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
    
case OPTIC_FLOW_SPEED 
    OpticFlowSpeed_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, skip_trials);
   
case LOOMING_SPEED
    LoomingSpeed_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, skip_trials);
    
case OPTIC_FLOW_COHER 
    OpticFlowCoherence_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE, skip_trials);

case OPTIC_FLOW_RGBLUM 
    OpticFlowContrast_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin,  StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);
   
case OPTIC_FLOW_EL
    OpticFlowElavation_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);

case OPTIC_FLOW_AZ    
    OpticFlowAzimuth_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);

case OPTIC_FLOW_GASSIAN_AMP
    OpticFlowGassianAmp_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);

case LOOMING_GASSIAN_AMP
    LoomingGassianAmp_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffsetBin, StopOffsetBin, StartEventBin, StopEventBin, PATH, FILE);

%-------------------------------------------
% Fixed Part Ended
%-------------------------------------------
    
otherwise
    set(ListHandle, 'String','The Protocol # is not known.  Check the list of Protocol codes in TEMPO_Defs.m');
end

return;