%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SpikeSort2.m
%
%   @Author:    Christopher Broussard
%   @Date:      November, 2001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = SpikeSort2(varargin)

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    
    % Initialization
    handles.ProgramName = 'SpikeSort2';
    handles.isFileOpen = 0;     % Indicate no spike2 file is open.
    handles.plotThreshSpikes = 0;
    clear global spsData;

	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.




% --------------------------------------------------------------------
function varargout = FileMenu_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = ExitItem_Callback(h, eventdata, handles, varargin)
    close(gcf);
% --------------------------------------------------------------------
function varargout = ToolsMenu_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function varargout = BatchItem_Callback(h, eventdata, handles, varargin)
    spikesort_batch;




% -------------------------------------------------------------------------
%   Opens up a browse box that lets the user choose which data file
%   to load.
% -------------------------------------------------------------------------
function varargout = BrowseButton_Callback(h, eventdata, handles, varargin)
    dfltDir = 'Z:/Data';   % Default data directory
    dfltSuffix = '*.smr';      % Default file type.
    cd(dfltDir);
    
    % Locate the data file and check if anything was chosen.
    [dataFileName, dataPath] = uigetfile(dfltSuffix, 'Choose Data File');
    if (dataFileName == 0)
        return;
    end
    
    % Load the data path edit box with the path and filename.
    set(handles.DataPathEdit, 'String', [dataPath, dataFileName]);
    
    handles.dataFileName = dataFileName;
    
    % Determine if we're looking at a MOOG monkey.  MOOG monkeys have a
    % different directory structure so we have to handle that accordingly.
    if findstr(dataPath, 'Z:\Data\MOOG')
        handles.isMoogMonkey = 1;
    else
        handles.isMoogMonkey = 0;
    end
    
    % MOOG monkeys have a different directory structure than other monkeys.
    if handles.isMoogMonkey
        slashIndex = findstr(dataPath, '\');
        handles.monkeyName = dataPath(slashIndex(3)+1:slashIndex(4)-1);
    else
        handles.monkeyName = dataPath(13:length(dataPath) - 1);
    end
    
    set(handles.OpenPlotButton, 'Enable', 'on');
    
    guidata(h, handles);   % Make sure everything is stored.
    


% ---------------------------------------------------------------------------
% Close the main window procedure.
% ---------------------------------------------------------------------------
function varargout = figure1_CloseRequestFcn(h, eventdata, handles, varargin)
    % Make sure that the spike2 data file is closed.
    if (handles.isFileOpen == 1)
        cedfunction('SonCloseFile', handles.spike2file);
        disp('Spike2 data file closed');
        ChannelCleanup;     % Clean up wasted memory.
    end
    
    disp([handles.ProgramName, ' terminated']);
    delete(gcf);
    
    
% -------------------------------------------------------------------------
%   Cleans up any channel globals so memory won't be wasted.
% -------------------------------------------------------------------------
function ChannelCleanup()
    % Remove all the channel data from memory
    clear global CHAN1;
    clear global CHAN2;
    clear global CHAN32;
    clear global spsData;
    
    % Also clean up all the functions in memory.
    clear functions;
    
    disp('Channel data wiped');


% ---------------------------------------------------------------------------
%   Opens the .smr file and plots it.
% ---------------------------------------------------------------------------
function varargout = OpenPlotButton_Callback(h, eventdata, handles, varargin)
    handles.wdata = 0;
    
    % Make sure we don't plot any thresholded spikes yet.
    handles.plotThreshSpikes = 0;
    
    % This value keeps track of the # of exports we've done per file.
    handles.spikeSet = 0;
    
    % Get the file to open.
    fileName = get(handles.DataPathEdit, 'String');
    
    % Check to see if a file was given.
    if (isempty(fileName))
        errordlg('You must choose a file to open!', 'File Error');
        return;
    end
    
    % Close any old open files.
    if (handles.isFileOpen == 1)
        cedfunction('SonCloseFile', handles.spike2file);
        disp('Old file closed');
    end
    
    % Look at the corresponding .htb file and load the pre/post buffer times.
    htdFile = handles.dataFileName;
    htdFile(length(htdFile) - 2 : length(htdFile)) = 'htb';
    if handles.isMoogMonkey
        handles.htdFile = ['Z:/Data/Moog/', handles.monkeyName, '/Raw/', htdFile];
    else
        handles.htdFile = ['Z:/Data/Tempo/', handles.monkeyName, '/Raw/', htdFile];
    end
    fid = htbOpen(handles.htdFile);     % Open the htb file.
    ndbs = htbCount(fid);               % Find out the # of databases it holds.
    % Find the Events database and get the times.
    for (i = 1:ndbs)
        hd = htbGetHd(fid, i);   % Get the database header.
        
        if (strcmp(hd.title, 'Events'))
            hertz = hd.speed_units / hd.speed;    % see p366 in TEMPO v9 manual
   	        binwidth = (hd.skip + 1) / hertz;
	        epoch_start = hd.offset * binwidth;
   	        epoch_period = hd.period * binwidth;
            set(handles.PreEventBufferEdit, 'String', num2str(epoch_start));
            set(handles.PostEventBufferEdit, 'String', num2str(-epoch_start + epoch_period));
            handles.PreEventBuffer = epoch_start;
            handles.PostEventBuffer = -epoch_start + epoch_period;
            break;
        end
    end
    htbClose(fid);
    
    % Open up the Spike2 data file.
    handles.spike2file = cedfunction('SonOpenOldFile', fileName, 1);
    handles.isFileOpen = 1;
    disp('File opened successfully');
    
    % Load ADC channel 1 into memory.
    global CHAN1;
    set(handles.ProcessText, 'String', 'Processing Channel 1');   % Update the status label
    handles.maxTime1 = cedfunction('SonChanMaxTime', handles.spike2file, 0);
    handles.chand = cedfunction('SonChanDivide', handles.spike2file, 0);
    numpts = cedfunction('SonGetNumADCPts', handles.spike2file, 0, 10000, 0, handles.maxTime1, handles.chand);
    [CHAN1, num] = cedfunction('SonGetAllADCData', handles.spike2file, 0, numpts, 0, handles.maxTime1, handles.chand);
    
    % Load ADC Channel 2 into memory (Window Discriminator).
    global CHAN2;
    if (get(handles.LoadDiscrimSpikesCheck, 'Value'))
        set(handles.ProcessText, 'String', 'Processing Channel 2');   % Update the status label
        handles.maxTime2 = cedfunction('SonChanMaxTime', handles.spike2file, 1);
        if (handles.maxTime2 > 0)
            handles.wdata = 1;
            [num, CHAN2] = cedfunction('SonGetEventData', handles.spike2file, 1, 1000000, 0, handles.maxTime2);
        else
            disp('No Channel 2');
        end
    else
        handles.maxTime2 = 0;
    end
    
    % Load Event channel 32 into memory.
    global CHAN32;
    set(handles.ProcessText, 'String', 'Processing Channel 32');  % Update the status label
    handles.maxTime32 = cedfunction('SonChanMaxTime', handles.spike2file, 31);
    [num, CHAN32] = cedfunction('SonGetMarkData', handles.spike2file, 31, 1000000, 0, handles.maxTime32);
    
    % Fix the Event channel data if necessary.  If the file is unfixable,
    % throw up an error message and close the file.
    [fixable, errsFound] = FixSMRData(handles);
    if (fixable == 0)
        errordlg('Could not fix SMR file, closing file!', 'SMR Fix Error');
        ClosePlotButton_Callback(h, eventdata, handles, varargin);
        return;
    end
    
    % Determine max and min range values for the voltages.
    maxVoltage = abs(double(max(CHAN1)) / (2^15 - 1) * 5);
    minVoltage = abs(double(min(CHAN1)) / (2^15 - 1) * 5);
    handles.rangeMax = round(max([maxVoltage, minVoltage])) + .5;
    
    % Set the upper and lower threshold sliders' min and max.
    set(handles.UpperboundSlider, 'Enable', 'on');
    set(handles.LowerboundSlider, 'Enable', 'on');
    set(handles.PlotSlider, 'Enable', 'on');
    set(handles.ApplyZoomButton, 'Enable', 'on');
    set(handles.UpperboundSlider, 'Min', 0);
    set(handles.UpperboundSlider, 'Max', handles.rangeMax);
    set(handles.LowerboundSlider, 'Min', -handles.rangeMax);
    set(handles.LowerboundSlider, 'Max', 0);
    set(handles.UpperboundSlider, 'Value', 0);
    set(handles.LowerboundSlider, 'Value', 0);
    get(handles.LowerboundSlider, 'Max');
    get(handles.UpperboundSlider, 'Max');
    
    % Set the slider step.
    set(handles.PlotSlider, 'SliderStep', [1 / (10 * GetCurrentZoomLevel(handles)), 1 / (2 * GetCurrentZoomLevel(handles))]);
    
    % Plot the data
    set(handles.ProcessText, 'String', 'Plotting...');  % Update the status label
    handles.currentZoom = GetCurrentZoomLevel(handles);
    [handles.startPoint, handles.endPoint, handles.upperLine, handles.lowerLine] = PlotSpikeData(handles);
    set(handles.upperLine, 'EraseMode', 'xor');
    set(handles.lowerLine, 'EraseMode', 'xor');
    
    % Update the Status string.
    set(handles.ProcessText, 'String', 'Finished Processing');
    
    % Disable the plot button and enable the close button.
    set(handles.OpenPlotButton, 'Enable', 'off');
    set(handles.ClosePlotButton, 'Enable', 'on');
    set(handles.AutoThreshButton, 'Enable', 'on');
    set(handles.LoadDiscrimSpikesCheck, 'Enable', 'off');
    set(handles.SortSpikesButton, 'Enable', 'on');
    
    % Set the sorted spikes output file.
    if handles.isMoogMonkey
       set(handles.OutputEdit, 'String', ['Z:\Data\Moog\', handles.monkeyName, '\CED\SortedSpikes\', handles.dataFileName(1:length(handles.dataFileName) - 3), 'mat']); 
    else
        set(handles.OutputEdit, 'String', ['Z:\Data\CED\', handles.monkeyName, '\SortedSpikes\', handles.dataFileName(1:length(handles.dataFileName) - 3), 'mat']);
    end
    
    % Calculate the single unit spontaneous level from the window
    % discriminator data.
    if handles.maxTime2 > 0
        handles.WinDiscrimSpontRate = FindWinDiscrimSpontRate(handles);
        set(handles.SpontRateLabel, 'String', num2str(handles.WinDiscrimSpontRate));
    else
        handles.WinDiscrimSpontRate = 0;
    end
    
    guidata(h, handles);
    
    
    
% --------------------------------------------------------------------
function currentZoom = GetCurrentZoomLevel(handles)
    % Determine the zoom level.
    str = get(handles.AutoZoomPopup, 'String');
    index = get(handles.AutoZoomPopup, 'Value');
    currentZoom = str2num(str{index});
    
    return;


% ----------------------------------------------------------------------------
%   Applies the zoom level and replots the data.
% ----------------------------------------------------------------------------
function varargout = ApplyZoomButton_Callback(h, eventdata, handles, varargin)
    
    handles.currentZoom = GetCurrentZoomLevel(handles);
    % Set the slider step.
    set(handles.PlotSlider, 'SliderStep', [1 / (10 * GetCurrentZoomLevel(handles)), 1 / (2 * GetCurrentZoomLevel(handles))]);
    [handles.startPoint, handles.endPoint, handles.upperLine, handles.lowerLine] = PlotSpikeData(handles);
    set(handles.upperLine, 'EraseMode', 'xor');
    set(handles.lowerLine, 'EraseMode', 'xor');
    
    guidata(h, handles);


% --------------------------------------------------------------------
function varargout = ClosePlotButton_Callback(h, eventdata, handles, varargin)
    % Close the .smr file.
    cedfunction('SonCloseFile', handles.spike2file);
    disp('Spike2 data file closed');
    
    % Clear the axes and the channel data.
    cla;
    clear global CHAN1;
    clear global CHAN2;
    clear global CHAN32;
    clear global spsData;
    disp('Axes cleared and channel data wiped');
    handles.isFileOpen = 0;
    
    set(h, 'Enable', 'off');
    set(handles.ApplyZoomButton, 'Enable', 'off');
    set(handles.ExportButton, 'Enable', 'off');
    set(handles.OpenPlotButton, 'Enable', 'on');
    set(handles.PlotSlider, 'Value', 0.0);
    set(handles.DataPathEdit, 'String', '');
    set(handles.ProcessText, 'String', '');
    set(handles.UpperboundSlider, 'Enable', 'off');
    set(handles.LowerboundSlider, 'Enable', 'off');
    set(handles.PlotSlider, 'Enable', 'off');
    set(handles.UpperBoundText, 'String', '');
    set(handles.LowerBoundText, 'String', '');
    set(handles.OutputEdit, 'String', '');
    set(handles.LoadDiscrimSpikesCheck, 'Enable', 'on');
    set(handles.SortSpikesButton, 'Enable', 'off');
    set(handles.SortList, 'Value', 1);
    set(handles.SortList, 'String', '');
    set(handles.AutoThreshButton, 'Enable', 'off');
    set(handles.OpenPlotButton, 'Enable', 'off');
    
    guidata(h, handles);


% -----------------------------------------------------------------------
% Adjust the viewing area of the plot taking into account the current
% zoom level.
% -----------------------------------------------------------------------
function varargout = PlotSlider_Callback(h, eventdata, handles, varargin)
    [handles.startPoint, handles.endPoint, handles.upperLine, handles.lowerLine] = PlotSpikeData(handles);
    set(handles.upperLine, 'EraseMode', 'xor');
    set(handles.lowerLine, 'EraseMode', 'xor');
    
    guidata(h, handles);



% -----------------------------------------------------------------------------
%   Sets the upper threshold limit.
% -----------------------------------------------------------------------------
function varargout = UpperboundSlider_Callback(h, eventdata, handles, varargin)
    ubound = get(handles.UpperboundSlider, 'Value');
    set(handles.upperLine, 'XData', [handles.startPoint, handles.endPoint], 'YData', [ubound, ubound]);
    
    % Update the text field that displays the current value of the slider.
    set(handles.UpperBoundText, 'String', num2str(get(h, 'Value')));
    
    guidata(h, handles);


% -----------------------------------------------------------------------------
%   Sets the lower threshold limit.
% -----------------------------------------------------------------------------
function varargout = LowerboundSlider_Callback(h, eventdata, handles, varargin)
    lbound = get(handles.LowerboundSlider, 'Value');
    set(handles.lowerLine, 'XData', [handles.startPoint, handles.endPoint], 'YData', [lbound, lbound]);
    
    % Update the text field that displays the current value of the slider.
    set(handles.LowerBoundText, 'String', num2str(get(h, 'Value')));
    
    guidata(h, handles);



% -------------------------------------------------------------------------
%   Exports thresholded ADC and Event data to a .mat file.  All imported
%   data is stuffed into an array of structures called 'spsData'.
% -------------------------------------------------------------------------
function varargout = ExportButton_Callback(h, eventdata, handles, varargin)
    global spsData;
       
    % Dump "spsData" to a .mat file.
    fileName = get(handles.OutputEdit, 'String');
    wdata = handles.wdata;
    eval(['save ', fileName, ' spsData wdata']);
    
    set(handles.ProcessText, 'String', 'Export Complete');
    
    guidata(h, handles);
    


% --------------------------------------------------------------------
function varargout = SortSpikesButton_Callback(h, eventdata, handles, varargin)
    % Set the values for the upper and lower thresholds.
    handles.utValue = get(handles.UpperboundSlider, 'Value');
    handles.ltValue = get(handles.LowerboundSlider, 'Value');
    
    % Sort the spikes.
    handles = SortSpikes(handles, 0);
    
    guidata(h, handles);


    
% -----------------------------------------------------------------------------
%   Automatically sets the threshold based on the average amount of
%   spontaneity.
% -----------------------------------------------------------------------------
function varargout = AutoThreshButton_Callback(h, eventdata, handles, varargin)
    codeTimes = FindSpontaneousCodeTimes(handles);
    
    threshLevel = str2num(get(handles.SpontEdit, 'String'));
    
    % Add in the window discriminator sontaneous rate if checked.
    if get(handles.AddSUSpontRateBox, 'Value')
        threshLevel = threshLevel + handles.WinDiscrimSpontRate;
    end
    
    % Find upper threshold level.
    if (get(handles.UpperThreshCheck, 'Value'))
        tValue = FindThreshold(codeTimes, 2, handles, handles.rangeMax / 2, threshLevel, 1, 0);
        set(handles.UpperboundSlider, 'Value', tValue);
        set(handles.upperLine, 'XData', [handles.startPoint, handles.endPoint], 'YData', [tValue, tValue]);
    
        % Update the text field that displays the current value of the slider.
        set(handles.UpperBoundText, 'String', num2str(get(handles.UpperboundSlider, 'Value')));
    end
    
    % Find lower threshold level.
    if (get(handles.LowerThreshCheck, 'Value'))
        tValue = FindThreshold(codeTimes, 2, handles, handles.rangeMax / 2, threshLevel, -1, 0);
        set(handles.LowerboundSlider, 'Value', -tValue);
        set(handles.lowerLine, 'XData', [handles.startPoint, handles.endPoint], 'YData', [-tValue, -tValue]);
    
        % Update the text field that displays the current value of the slider.
        set(handles.LowerBoundText, 'String', num2str(get(handles.LowerboundSlider, 'Value')));
    end
    
    set(handles.ProcessText, 'String', 'Threshold found');
    
    guidata(h, handles);
        
    
    
% --------------------------------------------------------------------------------
%   Clears all channel sorts out of memory and the Sort List.
% --------------------------------------------------------------------------------
function varargout = ClearSortListButton_Callback(h, eventdata, handles, varargin)
    clear global spsData;   % Trash all spsData.
    handles.spikeSet = 0;   % Reset which spike set we are on.
    set(handles.SortList, 'Value', 1);
    set(handles.SortList, 'String', '');
    set(handles.ExportButton, 'Enable', 'off');
    
    guidata(h, handles);
    
    
    
% --------------------------------------------------------------------------------
%   Removes a sort from spsData.
% --------------------------------------------------------------------------------
function varargout = RemoveSortButton_Callback(h, eventdata, handles, varargin)
    global spsData;
    
    % Find the index of the list item to be deleted.
    listIndex = get(handles.SortList, 'Value');
    listElements = get(handles.SortList, 'String');
    
    % Remove an element out of the spsData list and SortList.
    for (i = listIndex : length(spsData) - 1)
        spsData(i) = spsData(i + 1);
        listElements(i) = listElements(i + 1);
    end
    spsData(length(spsData)) = [];
    listElements(length(listElements)) = '';
    set(handles.SortList, 'String', listElements);
    handles.spikeSet = handles.spikeSet - 1;
    
    % Set which item in the SortList is highlighted.
    if (length(spsData) == 0)
        set(handles.SortList, 'Value', 1);
    elseif (length(spsData) < listIndex)
        set(handles.SortList, 'Value', length(spsData));
    else
        set(handles.SortList, 'Value', listIndex);
    end
    
    % Disable the remove and clear buttons if we've removed everything out of spsData.
    if (length(listElements) == 0)
        handles.spikeSet = 0;
        clear global spsData;
        set(handles.RemoveSortButton, 'Enable', 'off');
        set(handles.ClearSortListButton, 'Enable', 'off');
        set(handles.ExportButton, 'Enable', 'off');
    end
    
    guidata(h, handles);
    
    

% --------------------------------------------------------------------------
% Do nothing functions.  They're only defined here so Matlab will stop
% complaining.
% --------------------------------------------------------------------------
function varargout = AutoZoomPopup_Callback(h, eventdata, handles, varargin)
function varargout = UpperThreshCheck_Callback(h, eventdata, handles, varargin)
function varargout = LowerThreshCheck_Callback(h, eventdata, handles, varargin)
function varargout = OutputEdit_Callback(h, eventdata, handles, varargin)
function varargout = OutputBrowseButton_Callback(h, eventdata, handles, varargin)
function varargout = SpontEdit_Callback(h, eventdata, handles, varargin)
function varargout = LoadDiscrimSpikesCheck_Callback(h, eventdata, handles, varargin)
function varargout = SortList_Callback(h, eventdata, handles, varargin)


function varargout=AddSUSpontRateBox_Callback(h,eventdata,handles, varargin)


