%This function attempts to pool all five useful gratings analyses' data
%together for use in calculating latency RLS 1-20-08
function Pool_Data(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);   

global pooling pool_file pool_data1 pool_data2 pool_size pool_data pool_trials bin_size window_size, 
bin_size = .003;
window_size = .021;
pool_data1 = [];
pool_data2 = [];
pool_size = 0;
pool_trials = 0;
pooling = 1;
first_analysis = FILE(10);
if first_analysis == '2';
    SF_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, 1);
elseif first_analysis == '3';
    TF_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, 1);   
elseif first_analysis == '4';
    Orientation_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, 1);
elseif first_analysis == '5';
    Contrast_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, 1); 
elseif first_analysis == '6';
    Size_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, 1); 
end

for i = 2:6
    if int2str(i) ~= first_analysis;
        FILE(10) = int2str(i);
        pool_file = FILE;
        TEMPO_GUI_Switchyard('file open');
        TEMPO_GUI_Switchyard('load data');
        TEMPO_GUI_Switchyard('analyze');
        %if i == 2;
        %    SF_PSTH(pool_data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, 1);
        %elseif i == 3;
        %    TF_PSTH(pool_data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, 1);   
        %elseif i == 5;
        %    Contrast_PSTH(pool_data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, 1); 
        %elseif i == 6;
        %    Size_PSTH(pool_data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, 1); 
        %end
    end
end

pooling = 0;
pool_data1 = [];
pool_data2 = [];
pool_size = 0;
pool_trials = 0;
return
           
           