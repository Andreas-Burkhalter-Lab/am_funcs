function varargout = SpikeSort_Batch(varargin)
% SPIKESORT_BATCH Application M-file for SpikeSort_Batch.fig
%    FIG = SPIKESORT_BATCH launch SpikeSort_Batch GUI.
%    SPIKESORT_BATCH('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 01-Dec-2005 16:30:41

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
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
function varargout = BatchFileEdit_Callback(h, eventdata, handles, varargin)


% ----------------------------------------------------------------------------
% Loads the batchfile into memory.
% ----------------------------------------------------------------------------
function varargout = LoadBatchButton_Callback(h, eventdata, handles, varargin)
    dfltDir = 'Z:/Data/Tempo/Batch Files';   % Default data directory
    dfltSuffix = '*.m';      % Default file type.
    cd(dfltDir);

    % Locate the batch file and check if anything was chosen.
    [batchFileName, batchPath] = uigetfile(dfltSuffix, 'Choose Batch File');
    if (batchFileName == 0)
        return;
    end
    
    set(handles.BatchFileEdit, 'String', [batchPath, batchFileName]);
    
    % Use textread to grab the path and filename of all files to be processed.
    [handles.filePath, handles.fileName] = textread([batchPath, batchFileName], '%s%s%*[^\n]', 'commentstyle', 'matlab');
    
    % Store the monkey's name.
    handles.monkeyName = batchPath(13:length(batchPath) - 1);
    
    guidata(h, handles);



% --------------------------------------------------------------------
function varargout = StartButton_Callback(h, eventdata, handles, varargin)
    handles.wData = 0;
    
    % Make sure these variables are clear.
    clear global CHAN1;
    clear global CHAN2;
    clear global CHAN32;
    clear global spsData;
    
    % Make sure the user selects an upper and/or lower threshold.
    if (get(handles.UpperThreshCheck, 'Value') == 0 & get(handles.LowerThreshCheck, 'Value') == 0)
        errordlg('Please select a threshold!', 'Error');
        return;
    end
    
    % Clear the FileList.
    set(handles.FileList, 'Value', 1);
    set(handles.FileList, 'String', '');

    buff = fopen('Z:\LabTools\Matlab\TEMPO_Analysis\spikesort_batch.log', 'w');
    
    for (i = 1:length(handles.fileName))
        fileList = get(handles.FileList, 'String');
        listLength = length(fileList);
        handles.spikeSet = 0;
        
        % Generate .smr filpath and filename.
        x = find(handles.filePath{i} == '\');
        if isempty(strfind(handles.filePath{i}, 'MOOG'))
            handles.smrFilePath = ['Z:/Data/CED/', handles.filePath{i}(x(3)+1 : x(4) - 1), '/'];
        else
            handles.smrFilePath = ['Z:/Data/MOOG/', handles.filePath{i}(x(3)+1 : x(4) - 1), '/CED/'];
        end
        handles.smrFileName = [handles.fileName{i}(1:length(handles.fileName{i}) - 4), '.smr'];
        handles.dataFileName = handles.smrFileName;
        
        % Add the file being processed to the list.
        fileList{listLength + 1} = ['Analyzing -- ', handles.smrFilePath, handles.smrFileName, ' ........... '];
        set(handles.FileList, 'String', fileList);
        set(handles.FileList, 'Value', listLength + 1);
        
        % Call all the processing functions.
        disp('Loading spike data');
        [handles, fixable] = LoadSpikeData(handles, i);
        
        % Get the spontaneuos level we're searching for.
        handles.usLevel = str2num(get(handles.UpperSLevel, 'String'));
        handles.lsLevel = str2num(get(handles.LowerSLevel, 'String'));
            
        % Add in the window discriminator spontaneous rate if checked.
        if get(handles.AddSUSpontRateBox, 'Value')
            handles.usLevel = handles.usLevel + handles.WinDiscrimSpontRate;
            handles.lsLevel = handles.lsLevel + handles.WinDiscrimSpontRate;
        end

        % Go to the next file if there was a problem opening the file.
        if (handles.spike2file < 0 | fixable == 0)
            % Show that the file was bad or corrupt.
            if (fixable == 0)
                fileList{listLength + 1} = [fileList{listLength + 1}, 'Corrupt File'];
            else
                fileList{listLength + 1} = [fileList{listLength + 1}, 'Bad File'];
            end
            % Update the file list.
            set(handles.FileList, 'String', fileList);
            set(handles.FileList, 'Value', listLength + 1);

            % write out string to log file to keep record of results
            fprintf(buff, '%s \n', fileList{listLength + 1});
            continue;
        end
        
        disp('Finding thresholds');
        handles = FindThresholds(handles);
        disp('Sorting data');
        handles = SortSpikes(handles, 1);
        
        % Close the .smr file.
        % cedfunction('SonCloseFile', handles.spike2file);
        
        % Export spsData
        global spsData;
        wdata = handles.wData;
        eval(['save ', handles.smrFilePath, '/SortedSpikes/', handles.smrFileName(1:length(handles.smrFileName) - 3), 'mat', ' spsData wdata']);
        
        % Clear all data channels.
        clear global CHAN1;
        clear global CHAN32;
        clear global spsData;
        
        % Show that the file is done being processed.
        fileList{listLength + 1} = [fileList{listLength + 1}, 'Finished'];
        
        % Update the file list.
        set(handles.FileList, 'String', fileList);
        set(handles.FileList, 'Value', listLength + 1);
        
        % write out string to log file to keep record of results
        fprintf(buff, '%s \n', fileList{listLength + 1});
    end
    
    fclose(buff);

    % Indicate that everything is finished.
    fl = get(handles.FileList, 'String');
    fl{length(fl) + 1} = 'Done';
    set(handles.FileList, 'String', fl);
    set(handles.FileList, 'Value', length(fl));
    
    guidata(h, handles);
    
    
    
% ----------------------------------------------------------------------------------
%   Finds the treshold voltage values specified in the batch file.
% ----------------------------------------------------------------------------------
function handles = FindThresholds(handles)
    codeTimes = FindSpontaneousCodeTimes(handles);
    
    % Find upper threshold level.
    if (get(handles.UpperThreshCheck, 'Value'))
        handles.utValue = FindThreshold(codeTimes, 2, handles, handles.rangeMax / 2, handles.usLevel, 1, 1);
        disp('Upperthreshold found');
    end
    
    % Find lower threshold level.
    if (get(handles.LowerThreshCheck, 'Value'))
        handles.ltValue = -FindThreshold(codeTimes, 2, handles, handles.rangeMax / 2, handles.lsLevel, -1, 1);
        disp('Lowerthreshold found');
    end
    
    return;
    
    
    
% ----------------------------------------------------------------------------
%   Loads spike data from a file and processes it.
% ----------------------------------------------------------------------------
function [handles, fixable] = LoadSpikeData(handles, dumdex)
    fixable = 1;
    
    % Look at the corresponding .htb file and load the pre/post buffer times.
    handles.htdFile = [handles.filePath{dumdex}, handles.fileName{dumdex}];
    fid = htbOpen(handles.htdFile);     % Open the htb file.
    ndbs = htbCount(fid);    % Find out the # of databases it holds.
    % Find the Events database and get the times.
    for (i = 1:ndbs)
        hd = htbGetHd(fid, i);   % Get the database header.
        
        if (strcmp(hd.title, 'Events'))
            hertz = hd.speed_units / hd.speed;    % see p366 in TEMPO v9 manual
   	        binwidth = (hd.skip + 1) / hertz;
	        epoch_start = hd.offset * binwidth;
   	        epoch_period = hd.period * binwidth;
            handles.PreEventBuffer = epoch_start;
            handles.PostEventBuffer = -epoch_start + epoch_period;
            break;
        end
    end
    htbClose(fid);
    
    % Open up the Spike2 data file.
    handles.spike2file = cedfunction('SonOpenOldFile', [handles.smrFilePath, handles.smrFileName], 1);
    if (handles.spike2file >= 0)
        disp('File opened successfully');
    else
        disp('Could not open smr file.');
        return;
    end
    
    % Load ADC channel 1 into memory.
    global CHAN1;
    handles.maxTime1 = cedfunction('SonChanMaxTime', handles.spike2file, 0);
    handles.chand = cedfunction('SonChanDivide', handles.spike2file, 0);
    numpts = cedfunction('SonGetNumADCPts', handles.spike2file, 0, 10000, 0, handles.maxTime1, handles.chand);
    [CHAN1, num] = cedfunction('SonGetAllADCData', handles.spike2file, 0, numpts, 0, handles.maxTime1, handles.chand);
    
    % Load ADC Channel 2 into memory (Window Discriminator).
    global CHAN2;
    handles.maxTime2 = cedfunction('SonChanMaxTime', handles.spike2file, 1);
    if (handles.maxTime2 > 0)
        [num, CHAN2] = cedfunction('SonGetEventData', handles.spike2file, 1, 1000000, 0, handles.maxTime2);
        handles.wData = 1;
    else
        disp('No Channel 2');
        handles.maxTime2 = 0;
        handles.wData = 0;
    end
    
    % Load Event channel 32 into memory.
    global CHAN32;
    handles.maxTime32 = cedfunction('SonChanMaxTime', handles.spike2file, 31);
    [num, CHAN32] = cedfunction('SonGetMarkData', handles.spike2file, 31, 1000000, 0, handles.maxTime32);
    
    % Determine max and min range values for the voltages.
    maxVoltage = abs(double(max(CHAN1)) / (2^15 - 1) * 5);
    minVoltage = abs(double(min(CHAN1)) / (2^15 - 1) * 5);
    handles.rangeMax = round(max([maxVoltage, minVoltage])) + .5;
    
    % Close the open spike file.
    cedfunction('SonCloseFile', handles.spike2file);
    
    % We call this function to fix any errors in CHAN32.
    disp('Checking Data Integrity...');
    [fixable, errsFound] = FixSMRData(handles);
    
    % Calculate the single unit spontaneous level from the window
    % discriminator data.
    if handles.maxTime2 > 0
        handles.WinDiscrimSpontRate = FindWinDiscrimSpontRate(handles);
    else
        handles.WinDiscrimSpontRate = 0;
    end
    
    disp(['SU spont rate: ', num2str(handles.WinDiscrimSpontRate)]);
    
    return;
    


% --------------------------------------------------------------------
function varargout = UpperThreshCheck_Callback(h, eventdata, handles, varargin)
function varargout = LowerThreshCheck_Callback(h, eventdata, handles, varargin)


% --- Executes on button press in AddSUSpontRateBox.
function AddSUSpontRateBox_Callback(hObject, eventdata, handles)
% hObject    handle to AddSUSpontRateBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AddSUSpontRateBox


