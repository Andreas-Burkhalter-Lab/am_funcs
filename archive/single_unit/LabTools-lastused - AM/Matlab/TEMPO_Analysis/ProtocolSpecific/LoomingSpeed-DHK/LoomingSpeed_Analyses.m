%%%% 8/29/16 AM added optional 'skip_trials' argin for 'Plot Tuning Curve' - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped
%%%% 'skip_trials' are relative to beginning and end of run, not BegTrial

function LoomingSpeed_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, skip_trials)

	switch(Analysis{1})
	   case 'Plot Tuning Curve'
         LoomingSpeedTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials);
      case 'Plot Histograms'
         LoomingSpeed_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
      case 'Plot Auto and Cross Correlograms'
         CrossCorr(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
     case 'Fit Gaussian'
         LoomingSpeedTuningCurveGaussian(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

    end

return;