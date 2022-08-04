%%%% 8/29/16 AM added optional 'skip_trials' argin for 'Fit 2D Gaussian' - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped
%%%% 'skip_trials' are relative to beginning and end of run, not BegTrial

function GratingRFMap_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, skip_trials)

	switch(Analysis{1})
   	case 'Fit 2D Gaussian'
 	     GratingRFProfile(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials);           
    case 'Plot PSTH'
         PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Plot Spike Rasters'
         PlotRasters(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, Protocol);
    end

return;