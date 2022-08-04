%%%% 8/29/16 AM added optional 'skip_trials' argin for 'Plot Tuning Curve' - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped
%%%% 'skip_trials' are relative to beginning and end of run, not BegTrial

%%% AM edited 5/19/16: changed name from 'SizeTuning_Analyses' to differentiate
%%%     from SizeTuning_Analyses.m function written for grating size tuning

function DotSizeTuning_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials)

	switch(Analysis{1})
	   case 'Plot Tuning Curve'
           % AM 5/19/16: changed below function call from 'SizeTuningCurve' to specify
           % dot size tuning analysis rather than grating size tuning
           % analysis
         DotSizeTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials);
      case 'Plot Rasters/Histograms'
   		disp('plot rasters here');   
      case 'Plot Vergence Data' 
         disp('plot vergence here');   
      case 'Plot Auto and Cross Correlograms'
         CrossCorr(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);    
   end
return;