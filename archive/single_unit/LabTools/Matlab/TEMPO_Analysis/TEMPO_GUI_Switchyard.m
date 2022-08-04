%----------------------------------------------------------------------------------------------------------
%-- TEMPO_GUI_Switchyard.m: This function contains the callback routines for all of the TEMPO_GUI
%--	interface objects that require callback actions.  I have done it this way to keep all of the callback
%-- 	code in one centralized location.  GCD, 1/3/2000
%----------------------------------------------------------------------------------------------------------
function TEMPO_GUI_Switchyard(action)

global good_data;
global bad_data;
global listtext;
global EyeFlag;
global pool_data1; 
global pool_data2; 
global pool_size;
global pooling pool_file

TEMPO_Defs;
Path_Defs;

switch(action)
case 'file open'
    if pooling == 0;
        PathHandle = findobj(gcbf, 'Tag', 'Pathname');
        PATH = get(PathHandle, 'String');
        t = sprintf('%s*.htb', PATH);
        [fname, pname] = uigetfile(t,'Select .htb file');	
        if (fname ~= 0)	
            set(PathHandle, 'String', pname);
            FileHandle = findobj(gcbf, 'Tag', 'Filename');
            set(FileHandle, 'String', fname);
        end
    else
        PathHandle = findobj(gcbf, 'Tag', 'Pathname');
        PATH = get(PathHandle, 'String');
        FileHandle = findobj(gcbf, 'Tag', 'Filename');
        set(FileHandle, 'String', pool_file);
    end
case 'file info'
    PathHandle = findobj(gcbf, 'Tag', 'Pathname');
    PATH = get(PathHandle, 'String');
    FileHandle = findobj(gcbf, 'Tag', 'Filename');
    FILE = get(FileHandle, 'String');
    [return_val, info_text] = GetHTBInfo(PATH,FILE);
    ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
    if (return_val == -1)  %problem opening file
        listtext{length(listtext)+1} = 'ERROR: File could not be opened';
        set(ListHandle, 'String', listtext);
    else
        set(ListHandle, 'String', info_text);
    end
case 'load data'
    PathHandle = findobj(gcbf, 'Tag', 'Pathname');
    PATH = get(PathHandle, 'String');
    FileHandle = findobj(gcbf, 'Tag', 'Filename');
    FILE = get(FileHandle, 'String');
    ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
    listtext = [];
    listtext{1} = 'Loading Data...';
    set(ListHandle, 'String', listtext);
    [return_val, good_data, bad_data, newChans] = LoadTEMPOData(PATH,FILE);
    %% note that newChans ~= 0 if Packaroni loads new spike channels - BJP
    
    if (return_val == -1)  %problem opening file
        listtext{length(listtext)+1} = 'ERROR: File could not be opened';
        set(ListHandle, 'String', listtext);
    else	%data are loaded, now setup some interface items
        listtext{1} = 'Loading Data...DONE';
        set(ListHandle, 'String', listtext);
        ProtocolHandle = findobj(gcbf, 'Tag', 'Protocol');
        str = protocol_names{good_data.one_time_params(PROTOCOL)+1};		%get the protocol string description
        set(ProtocolHandle, 'String', str);		% show the protocol string description
        set(ProtocolHandle, 'Value', good_data.one_time_params(PROTOCOL));
        AnalysisHandle = findobj(gcbf, 'Tag', 'AnalysisPopup');
        set(AnalysisHandle, 'String', analysis_strings{good_data.one_time_params(PROTOCOL)+1});
        %BJP 1/14/01 - I included the next line so that analysis options not disappear if going from many to few analyses
        %   Otherwise, if going from lots of options to few, analysis drop menue disappears
        set(AnalysisHandle, 'Value', 1);
        if (pooling == 1 & FILE(10) == '2')
            set(AnalysisHandle, 'Value', 6);
        elseif pooling ==1
            set(AnalysisHandle, 'Value', 4);
        end
        %now, set up the Start Time and Stop Time Popup Menus
        StartHandle = findobj(gcbf, 'Tag', 'StartTimePopup');
        set(StartHandle, 'String', event_names);		
        set(StartHandle, 'Value', 4);		%default to Visual Stim ON
        StopHandle = findobj(gcbf, 'Tag', 'StopTimePopup');
        set(StopHandle, 'String', event_names);	
        set(StopHandle, 'Value', 5);		%default to Visual Stim OFF
        %set beginning and ending trial defaults        
        BegTrialHandle = findobj(gcbf, 'Tag', 'StartTrial');
        set(BegTrialHandle, 'String', '1');		
        EndTrialHandle = findobj(gcbf, 'Tag', 'StopTrial');
        set(EndTrialHandle, 'String', num2str(size(good_data.event_data,3)));		
        %set start and stop timing offsets for analysis - BJP 3/1/00
        StartOffsetHandle = findobj(gcbf, 'Tag', 'StartOffset');		
        set(StartOffsetHandle, 'String', '0');
        StopOffsetHandle = findobj(gcbf, 'Tag', 'StopOffset');
        set(StopOffsetHandle, 'String', '0');          
        %set-up spike-channel popup menu
        good_chan_num = size(good_data.spike_data);
        nchan = good_chan_num(1); %good_data.htb_header{SPIKE_DB}.nchannels;
        popup_text = [];
        for i=1:nchan
            spks = good_data.spike_data(i,:,:);
            num_events = sum(sum(spks));	%count up all spikes in this spike channel
            global DescriptionsExist;
            
            if (DescriptionsExist)
                if (isempty(good_data.desc{i}))
                    temp = sprintf('Recorded Spikes (Channel %d):  (%d events)', i, num_events);
                else
                    temp = [good_data.desc{i}, sprintf(', (%d events)', num_events)];
                end
            else
                %good_data.desc{i} = sprintf('Recorded Spikes (Channel %d):   (%d events)', i, num_events);
                %temp = [good_data.desc{i}];
                temp = sprintf('Recorded Spikes (Channel %d):  (%d events)', i, num_events);
            end
            
            popup_text{length(popup_text)+1} = temp;
        end
        % add channels from LFP database if present
        if ~isempty(good_data.lfp_data)
            nchan2 = size(good_data.lfp_data,1);
            for i=1:nchan2
                num_events = size(good_data.lfp_data,2);	%count up all spikes in this spike channel
                temp = sprintf('Recorded LFP (Channel %d):   (%d events)', i, num_events);
                popup_text{length(popup_text)+1} = temp;
            end            
        end    
        SpikeHandle = findobj(gcbf, 'Tag', 'SpikePopup');
        set(SpikeHandle, 'String', popup_text);
        SpikeHandle2 = findobj(gcbf, 'Tag', 'SpikePopup2');
        set(SpikeHandle2, 'String', popup_text);
        
        popup_text = [];
        % add eye channels from eye database if present
        if ~isempty(good_data.eye_data)
            nchan = size(good_data.eye_data,1);
            for i=1:nchan
                switch i
                case LEYE_H
                    temp = sprintf('Left Hor (%d)', i);
                case LEYE_V
                    temp = sprintf('Left Vert (%d)', i);
                case REYE_H
                    temp = sprintf('Right Hor (%d)', i);
                case REYE_V
                    temp = sprintf('Right Vert (%d)', i);
                 case DA_H
                    temp = sprintf('D/A Hor (%d)', i);
                case DA_V
                    temp = sprintf('D/A Vert (%d)', i);
               end
                popup_text{length(popup_text)+1} = temp;                
            end    
            temp = sprintf('Left Eye (%d-%d)', LEYE_H, LEYE_V);
            popup_text{length(popup_text)+1} = temp;
            temp = sprintf('Right Eye (%d-%d)', REYE_H, REYE_V);
            popup_text{length(popup_text)+1} = temp;
            temp = sprintf('Both Eyes');
            popup_text{length(popup_text)+1} = temp;
        end    
        EyeHandle = findobj(gcbf, 'Tag', 'EyePopup');
        set(EyeHandle, 'String', popup_text);        
    end
case 'pick good trials'
    BadHandle = findobj(gcbf, 'Tag', 'Bad Trials');
    set(BadHandle, 'Value', 0);
case 'pick bad trials'
    GoodHandle = findobj(gcbf, 'Tag', 'Good Trials');
    set(GoodHandle, 'Value', 0);
case 'use sync pulses'
    GoodHandle = findobj(gcbf, 'Tag', 'Sync Pulses');

    %set(GoodHandle, 'Value', 0);   
    
case 'analyze'
    %gather up some pieces of information to pass along to the analysis switchyard
    ProtocolHandle = findobj(gcbf, 'Tag', 'Protocol');
    Protocol = get(ProtocolHandle, 'Value');
    AnalysisHandle = findobj(gcbf, 'Tag', 'AnalysisPopup');
    Analysis = analysis_strings{get(ProtocolHandle, 'Value')+1}(get(AnalysisHandle, 'Value'));
    StartHandle = findobj(gcbf, 'Tag', 'StartTimePopup');
    StartCode = get(StartHandle, 'Value');
    StopHandle = findobj(gcbf, 'Tag', 'StopTimePopup');
    StopCode = get(StopHandle, 'Value');
    BegTrialHandle = findobj(gcbf, 'Tag', 'StartTrial');
    BegTrial = str2num(get(BegTrialHandle, 'String'));
    EndTrialHandle = findobj(gcbf, 'Tag', 'StopTrial');
    EndTrial = str2num(get(EndTrialHandle, 'String'));
    StartOffsetHandle = findobj(gcbf, 'Tag', 'StartOffset');			% next four lines added BJP 3/1/00
    StartOffset = str2num(get(StartOffsetHandle, 'String'));
    StopOffsetHandle = findobj(gcbf, 'Tag', 'StopOffset');
    StopOffset = str2num(get(StopOffsetHandle, 'String'));
    SyncPulseHandle = findobj(gcbf,'Tag', 'Sync Pulses');
    UseSyncPulses = get(SyncPulseHandle,'Value');   
    SpikeHandle = findobj(gcbf, 'Tag', 'SpikePopup');
    SpikeChan = get(SpikeHandle, 'Value');
    SpikeHandle2 = findobj(gcbf, 'Tag', 'SpikePopup2');
    SpikeChan2 = get(SpikeHandle2, 'Value');
    EyeFlagHandle = findobj(gcbf, 'Tag', 'EyePopup');
    %    EyeFlag = str2num(get(EyeFlagHandle, 'String'));
    EyeFlag = get(EyeFlagHandle, 'Value');

    PathHandle = findobj(gcbf, 'Tag', 'Pathname');
    PATH = get(PathHandle, 'String');
    FileHandle = findobj(gcbf, 'Tag', 'Filename');
    FILE = get(FileHandle, 'String');
    %depending on how the Good Trials/Bad Trials radio buttons are set, analyze the appropriate set of trials
    GoodHandle = findobj(gcbf, 'Tag', 'Good Trials');
    good_flag = get(GoodHandle, 'Value');
    if (good_flag)
        select_data = good_data;
    else
        select_data = bad_data;
    end
    
    select_data.eye_flag = EyeFlag;
    % pass everything to a switchyard function which will chain off to the apprpriate analysis routines
    %batch_flag tells analysis switchyard whether or not to process batch files - JDN 3/20/00
    batch_flag = 0;
    Analysis_Switchyard(select_data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag, UseSyncPulses);
case 'view codes'   
    %depending on how the Good Trials/Bad Trials radio buttons are set, select the appropriate set of trials
    GoodHandle = findobj(gcbf, 'Tag', 'Good Trials');
    good_flag = get(GoodHandle, 'Value');
    if (good_flag)
        select_data = good_data;
    else
        select_data = bad_data;
    end
    %for each selected trial, find all of the non-zero elements of event_data; these are the code values
    %then sprintf them to a string and display them in the console window.
    num_trials = size(select_data.event_data, 3);
    for (i = 1:num_trials)
        tt = find( select_data.event_data(1,:,i) ~= 0 );
        str1 = sprintf('%d ', select_data.event_data(1,tt,i));
        str2 = sprintf('Trial #%4d: ', i);
        tstr{i} = [str2 str1];
    end
    ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
    set(ListHandle, 'String', tstr);
    
case 'view conditions'   
        %gather up some pieces of information to pass along to the analysis switchyard
    ProtocolHandle = findobj(gcbf, 'Tag', 'Protocol');
    Protocol = get(ProtocolHandle, 'Value');
    BegTrialHandle = findobj(gcbf, 'Tag', 'StartTrial');
    BegTrial = str2num(get(BegTrialHandle, 'String'));
    EndTrialHandle = findobj(gcbf, 'Tag', 'StopTrial');
    EndTrial = str2num(get(EndTrialHandle, 'String'));
    PathHandle = findobj(gcbf, 'Tag', 'Pathname');
    PATH = get(PathHandle, 'String');
    FileHandle = findobj(gcbf, 'Tag', 'Filename');
    FILE = get(FileHandle, 'String');

    %depending on how the Good Trials/Bad Trials radio buttons are set, select the appropriate set of trials
    GoodHandle = findobj(gcbf, 'Tag', 'Good Trials');
    good_flag = get(GoodHandle, 'Value');
    if (good_flag)
        select_data = good_data;
    else
        select_data = bad_data;
    end
    [conditions, UniqueConds, param, num_conditions] = RegenerateConditionList(select_data, BegTrial, EndTrial, PATH, FILE, Protocol);
    
    textline = 1;
    tstr{textline} = '**** Condition List ****';
    textline = textline + 1;
    tstr{textline} = 'Condition ';
    for modality = 1:length(param);
       tstr{textline} = [tstr{textline} '   ' param{modality} ];      
    end

    textline = textline + 1;    
    for cond = 1:num_conditions
        tstr{textline} = [sprintf('%2d',UniqueConds(cond,end)) '    '];
        tstr{textline} = [tstr{textline} sprintf('%6.3f     ', UniqueConds(cond,1:end-1)) ];    
        textline = textline + 1;  
    end

    tstr{textline} = ' ';
    textline = textline + 1;
    tstr{textline} = '**** Listed by Trial ****';
    textline = textline + 1;

    
    num_trials = size(select_data.event_data, 3);
    for (i = 1:num_trials)
        cond = conditions(end,i);
        str1 = sprintf('%6.3f ', UniqueConds(cond + 1,1:end - 1 )    );
        str2 = sprintf('Trial #%4d: %3d / ', i, cond);
        tstr{textline + i} = [str2 str1];
    end

    ListHandle = findobj(gcbf, 'Tag', 'MessageBox');
    set(ListHandle, 'String', tstr);
    
otherwise
    disp('unknown action');
end