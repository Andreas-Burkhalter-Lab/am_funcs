%%%% 8/29/16 AM added optional 'skip_trials' argin for 'Plot Tuning Curve' - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped
%%%% 'skip_trials' are relative to beginning and end of run, not BegTrial

function OpticFlowCoherence_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, skip_trials)
switch(Analysis{1})
    case 'Plot Tuning Curve'
        OptFlMotionCoherenceCurves(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials);
    case 'Plot PSTH'
        OptFlCoherence_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Spike Rasters'
        OptFlCoherenceRasters(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
end
return;