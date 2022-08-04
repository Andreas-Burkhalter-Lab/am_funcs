function PlotEyeTraces(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol)

TEMPO_defs;
h = data.htb_header{EYE_DB};	%for convenience
eye_bin_width = 1000*(h.skip + 1) / (h.speed_units / h.speed);


[msac, emt] = CalcMsac(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);    
[msac_st, msac_trl]=find(msac'==1);	% y = trl; x = index into emt column
[msac_end, msac_trl2]=find(msac'==-1);	% y = trl; x = index into emt column

emh= (squeeze(data.eye_data(REYE_H,:,BegTrial:EndTrial) ) )';
emv = (squeeze(data.eye_data(REYE_V,:,BegTrial:EndTrial) ) )';

%calculate mean eye position during vstim period - approx fixation pt location
fixpt_x = mean( mean(emh(:, floor(find(data.event_data(1,:,1) == VSTIM_ON_CD)/eye_bin_width )  :floor(find(data.event_data(1,:,1) == VSTIM_OFF_CD)/eye_bin_width)) ) );
fixpt_y = mean( mean(emv(:, floor(find(data.event_data(1,:,1) == VSTIM_ON_CD)/eye_bin_width )  :floor(find(data.event_data(1,:,1) == VSTIM_OFF_CD)/eye_bin_width)) ) );

%offset eye positions relative to fixation point location
emh = emh - fixpt_x;
emv = emv - fixpt_y;

for trial = 1:size(emh,1)
    emh(trial, 1:floor(StartEventBin(trial)/eye_bin_width)  )= 0;
    emv(trial, 1:floor(StartEventBin(trial)/eye_bin_width)  ) = 0;
    emh(trial, (floor(StopEventBin(trial)/eye_bin_width + 1): end   )  ) = 0;
    emv(trial, (floor(StopEventBin(trial)/eye_bin_width + 1): end   )  ) = 0;
end



%emt = ( 1:eye_bin_width: length(emh) * eye_bin_width ) - find(data.event_data(1,:,1) == StartCode);

plot_flag = 0;
num_columns = 2;
num_rows = 5;
% plotting loops for traces and microsaccades
    for page =1:ceil( (EndTrial-BegTrial+1)/10)
        if (plot_flag == 1)
         figure
         subplot(num_rows, num_columns,1);
        end
        for col = 1:num_columns
            for row = 1:num_rows
                trial = (page - 1)*10 + row + (col - 1)*num_rows;
                msac_markers = [msac_st(find(msac_trl == trial))' msac_end(find(msac_trl2 == trial))']';  
                if (trial <= size(emh,1) )
                    %align time bins based sync pulse start signal
                    emt = ( 1:eye_bin_width: length(emh) * eye_bin_width ) - StartEventBin(trial);
                    if (plot_flag == 1)

                        subplot(num_rows, num_columns, (row*2 - 1*(2 - col) ) );
                        plot(emt, emh(trial,:), 'r-', emt, emv(trial,:),'b-'  );  
                        hold on;
                        plot([emt(msac_markers); emt(msac_markers) ], [-0.5*ones(length(msac_markers), 1) -1*ones(length(msac_markers), 1) ]', 'k-' );
                        ylim([-2 2]);
                        xlim([ (find(data.event_data(1,:,1) == StartCode)-find(data.event_data(1,:,1) == StartCode)) (find(data.event_data(1,:,1) == StopCode)- find(data.event_data(1,:,1) == StartCode) ) ]);
                        hold on;
                        text(-500,0,['Trial ' num2str(BegTrial + trial - 1)]);
                    end    
                end        
                msac_start_times{trial} = emt(msac_st(find(msac_trl == trial))');
                msac_end_times{trial} = emt(msac_end(find(msac_trl2 == trial))' );
            end
        end
    end

output = 1;

    %output tuning curve metrics
    if (output == 1)
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end   
        PATHOUT = [PATH(1:i) 'Analysis\Saccades\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        FILEOUT = [FILE(1:i) 'msc'];
        
        eval(['save ' PATHOUT FILEOUT  ' msac_start_times msac_end_times StartCode StopCode BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
    end
