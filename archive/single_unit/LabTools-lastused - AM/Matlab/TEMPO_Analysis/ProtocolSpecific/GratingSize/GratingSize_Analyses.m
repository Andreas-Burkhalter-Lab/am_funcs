%%%% 8/29/16 AM added optional 'skip_trials' argin for 'Plot Tuning Curve' - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped
%%%% 'skip_trials' are relative to beginning and end of run, not BegTrial

function GratingSize_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, skip_trials)

global pool_size

	switch(Analysis{1})
   	case 'Plot Tuning Curve'
 	     SizeTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials);           
    case 'Plot PSTH'
         Size_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, 0);
    case 'Plot Spike Rasters'
         PlotSizeRasters(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
    case 'Pooling'
         if pool_size ~= 0;
             Size_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, 1);
         else
             Pool_Data(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
         end
    end

return;