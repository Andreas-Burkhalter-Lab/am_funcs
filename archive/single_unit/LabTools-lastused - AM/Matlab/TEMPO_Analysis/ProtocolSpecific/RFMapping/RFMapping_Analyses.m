%%%%%% AM note: rf mapping with drifting dots, not gratings
%%%% 8/29/16 AM added optional 'skip_trials' argin for 'Fith with 2D Gaussian' - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped
%%%% 'skip_trials' are relative to beginning and end of run, not BegTrial

function RFMapping_Analyses(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials)

	switch(Analysis{1})
        case 'Fit with 2D Gaussian'
            Gauss2D_Protocol(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials);
    end

return;