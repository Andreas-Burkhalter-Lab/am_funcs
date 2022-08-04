function BATCH_GUI_Tempo_Analysis(batchfiledir, filename, close_flag, print_flag, protocol_filter, analysis_filter)

TEMPO_Defs;
Path_Defs;

filename = [batchfiledir filename];
fid = fopen(filename);

line = fgetl(fid);
while (line ~= -1)
    
    %pause
    % format for batch files
    % PATH  FILE    
    
    if (line(1) ~= '%')

        % first remove any comment text at the end of a line (following a %), GCD, added 9/25/01
        comment_start = find(line == '%');
        if ~isempty(comment_start)
            line = line(1:(comment_start(1)-1));
        end

        spaces = isspace(line);
        space_index = find(spaces);
        
        %get path / file
        PATH = line(1:space_index(1) - 1);
        FILE = line(space_index(1) + 1:space_index(2) - 1)
        
        if protocol_filter ~= -1
            
            l = length(FILE);
            if (FILE(l-3:l) == '.htb')	% .htb extension already there
                filename = [PATH FILE];   %the HTB data file
                logfile = [PATH FILE(1:l-4) '.log'];   %the TEMPO log file
            else	%no extension in FILE, add extensions
                filename = [PATH FILE '.htb'];   %the HTB data file
                logfile = [PATH FILE '.log'];   %the TEMPO log file
            end
            
            %read in protocol type and add to cell array
            protocol_info = textread(logfile, '%s', 3);
            protocol_name = str2num(protocol_info{2});
            
            if protocol_name == protocol_filter
                %this protocol matches the ones that we want to analyze.
                %proceed with the analysis.
                do_file = 1;
            else %otherwise goto next file.
                do_file = 0;
            end
            
        else
            do_file = 1;
        end % if protocol_filter ~= 1 ...
        
        if do_file == 1
            %get analysis type
            if analysis_filter ~= -1
                TempAnalysis1 = analysis_strings{protocol_filter + 1}(analysis_filter);
                TempAnalysis = TempAnalysis1{1};
                i = 2;
                while(line(space_index(i)-1) ~= '''')
                    i = i+1;
                end
            else
                %special commands for analysis read
                i = 2;
                TempAnalysis = '';
                while(line(space_index(i)-1) ~= '''')
                    i = i+1;
                end
                
                TempAnalysis = strcat(TempAnalysis, line(space_index(2) + 2:space_index(i)-2));
            end
            
            Analysis = {};
            Analysis{1} = TempAnalysis;
            
            %get beginning trial
            BegTrial = line(space_index(i) + 1:space_index(i+1) - 1);
            BegTrial = str2num(BegTrial);
            
            %get ending trial
            EndTrial = line(space_index(i+1) + 1:space_index(i+2) - 1);
            EndTrial = str2num(EndTrial);
            
            %get startcode
            StartCode = line(space_index(i+2) + 1:space_index(i+3) - 1);
            StartCode = str2num(StartCode);
            
            %get startoffset
            StartOffset = line(space_index(i+3) + 1:space_index(i+4) - 1);
            StartOffset = str2num(StartOffset);
            
            %get stopcode
            StopCode = line(space_index(i+4) + 1:space_index(i+5) - 1);
            StopCode = str2num(StopCode);
            
            %get stopoffset
            StopOffset = line(space_index(i+5) + 1:space_index(i+6) - 1);
            StopOffset = str2num(StopOffset);
            
            %get select_data - good trials or bad trials
            good_flag = line(space_index(i+6) + 1:space_index(i+7) - 1);
            good_flag = str2num(good_flag);
            
            %get spikechannels & sync code flags
            
            if (i + 9 - 1 < length(space_index) )
                %if there are sync pulse value 
                SpikeChan = line(space_index(i+7) + 1:space_index(i+8) - 1);
                SpikeChan = str2num(SpikeChan);
                SpikeChan2 = line(space_index(i+8) + 1:space_index(i+9) - 1);
                SpikeChan2 = str2num(SpikeChan2);
                UseSyncPulses = line(space_index(i+9) +  1:length(line) - 1);
                UseSyncPulses = str2num(UseSyncPulses);        
            else
                SpikeChan = line(space_index(i+7) + 1:length(line));
                SpikeChan = str2num(SpikeChan);
                SpikeChan2 = SpikeChan;
                %if sync pulses flag not set, assume we want to use TEMPO START_CD & STOP_CD
                UseSyncPulses = 0;
            end   
            
            %load data file
            [return_val, g_data, b_data] = LoadTEMPOData(PATH,FILE);
            if (good_flag)
                select_data = g_data;
            else
                select_data = b_data;
            end
            
            if BegTrial == -1
                BegTrial = 1;		
            end
            
            if EndTrial == -1
                EndTrial = size(g_data.event_data,3);		
            end
            
            %get protocol type
            Protocol = g_data.one_time_params(PROTOCOL);		%get the protocol string descriptio
            
            %analyze data
            batch_flag = 1;
            if SpikeChan == 999
                chan_proc = size(select_data.spike_data);
                if chan_proc > 2
                    SpikeChan = chan_proc(1)
                    Analysis_Switchyard(select_data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag, UseSyncPulses);
                    if print_flag
                        printhandle = figure;
                        close(printhandle);
                        for printindex = 2:1:printhandle - 1
                            print(printindex);
                        end
                    end
                    if close_flag
                        closehandle = figure;
                        for closeindex = 2:1:closehandle
                            close(closeindex);
                        end
                    end
                else
                    t_string = sprintf('No MU Channels for %s.', FILE);
                    disp(t_string);
                end
            else
                Analysis_Switchyard(select_data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag, UseSyncPulses);
                if print_flag
                    printhandle = figure;
                    close(printhandle);
                    for printindex = 2:1:printhandle - 1
                        print(printindex);
                    end
                end
                if close_flag
                    closehandle = figure;
                    for closeindex = 2:1:closehandle
                        close(closeindex);
                    end
                end
            end
            
        end %do_file
    end %if (line(1) ~= '%')...
    
    line = fgetl(fid);
    
end %while
fclose(fid);