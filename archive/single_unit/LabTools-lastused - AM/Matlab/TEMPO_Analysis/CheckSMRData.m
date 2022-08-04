function varargout = CheckSMRData(varargin)
% CHECKSMRDATA M-file for CheckSMRData.fig
%      CHECKSMRDATA, by itself, creates a new CHECKSMRDATA or raises the existing
%      singleton*.
%
%      H = CHECKSMRDATA returns the handle to a new CHECKSMRDATA or the handle to
%      the existing singleton*.
%
%      CHECKSMRDATA('Property','Value',...) creates a new CHECKSMRDATA using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to CheckSMRData_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CHECKSMRDATA('CALLBACK') and CHECKSMRDATA('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CHECKSMRDATA.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CheckSMRData

% Last Modified by GUIDE v2.5 21-Aug-2002 14:17:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CheckSMRData_OpeningFcn, ...
    'gui_OutputFcn',  @CheckSMRData_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    varargout{1:nargout} = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CheckSMRData is made visible.
function CheckSMRData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for CheckSMRData
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CheckSMRData wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CheckSMRData_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function BatchFileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BatchFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function BatchFileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to BatchFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BatchFileEdit as text
%        str2double(get(hObject,'String')) returns contents of BatchFileEdit as a double


% -------------------------------------------------------------------------
% --- Executes on button press in StartButton.
% -------------------------------------------------------------------------
function StartButton_Callback(hObject, eventdata, handles)
% hObject    handle to StartButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global CHAN32;
totalFiles = 0;
totalCorrect = 0;
totalErrorsFound = 0;

% Clear the FileList.
set(handles.FileList, 'Value', 1);
set(handles.FileList, 'String', '');

% This iterates through each SMR file in the batch file, and checks it for
% errors.
for (i = 1:length(handles.fileName))
    fileList = get(handles.FileList, 'String');
    listLength = length(fileList);
    
    % Generate .smr filpath and filename.
    x = find(handles.filePath{i} == '\');
    handles.smrFilePath = ['Z:/Data/CED/', handles.filePath{i}(x(3)+1 : x(4) - 1), '/'];
    handles.smrFileName = [handles.fileName{i}(1:length(handles.fileName{i}) - 4), '.smr'];
    
    % Add the file being processed to the list.
    fileList{listLength + 1} = ['Job #', num2str(length(handles.fileName) - i + 1), ' -- ', handles.smrFilePath, handles.smrFileName, ' ............... '];
    set(handles.FileList, 'String', fileList);
    set(handles.FileList, 'Value', listLength + 1);
    
    % Call all the processing functions.
    handles = LoadSpikeData(handles, i);
    
    % Go to the next file if there was a problem opening the file.
    if (handles.spike2file < 0)
        % Show that the file was bad.
        fileList{listLength + 1} = [fileList{listLength + 1}, 'Bad File'];
        % Update the file list.
        set(handles.FileList, 'String', fileList);
        set(handles.FileList, 'Value', listLength + 1);
        continue;
    end
    
    % Derive the monkey name from the SMR file name.
    handles.monkeyName = handles.smrFilePath(13:length(handles.smrFilePath) - 1);
    
    % Make naming conventions jive with spikesort2.
    handles.dataFileName = handles.smrFileName;
    
    % Check the SMR data.
    totalFiles = totalFiles + 1;
    [fxbl, herrs] = FixSMRData(handles);
    totalCorrect = totalCorrect + fxbl;
    totalErrorsFound = totalErrorsFound + herrs;
    
    % Show that the file is done being processed.
    fileList{listLength + 1} = [fileList{listLength + 1}, 'Finished'];
    
    % Update the file list.
    set(handles.FileList, 'String', fileList);
    set(handles.FileList, 'Value', listLength + 1);
end % End -- for (i = 1:length(handles.fileName))

% Indicate that everything is finished.
fl = get(handles.FileList, 'String');
fl{length(fl) + 1} = 'Done';
set(handles.FileList, 'String', fl);
set(handles.FileList, 'Value', length(fl));

% Spit out some summary information.
fprintf('\n*************\nCheckSMRData Batch Summary:\n- No. Files Checked: %d\n- No. Trials w/ No Errors: %d\n', totalFiles, totalCorrect - totalErrorsFound);
fprintf('- No. Unrecoverable Files: %d\n- No. Corrected Files: %d\n', totalFiles - totalCorrect, totalErrorsFound);
fprintf('- Average Pct. Unrecoverable Files (bad files / total files): %.2f%%\n', (totalFiles - totalCorrect) / totalFiles * 100);
fprintf('- Average Pct. Files with Correctable Errors (file with fixable errors / total files): %.2f%%\n', totalErrorsFound / totalFiles * 100);

guidata(hObject, handles);


% --- Executes on button press in LoadButton.
function LoadButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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

fileList{1} = ['Number of files to process: ', num2str(length(handles.filePath))];
set(handles.FileList, 'String', fileList);
set(handles.FileList, 'Value', 1);
set(handles.StartButton, 'Enable', 'on');

guidata(hObject, handles);


% ----------------------------------------------------------------------------
%   Loads spike data from a file and processes it.
% ----------------------------------------------------------------------------
function handles = LoadSpikeData(handles, dumdex)    
% Look at the corresponding .htb file and load the pre/post buffer times.
handles.htdFile = [handles.filePath{dumdex}, handles.fileName{dumdex}];
fid = htbOpen(handles.htdFile);     % Open the htb file.
if (fid < 0)
    disp(['Could not open Tempo log: ', handles.htdFile]);
    handles.spike2file = -1;
    return;
else
    htbClose(fid);
end

% Open up the Spike2 data file.
handles.spike2file = cedfunction('SonOpenOldFile', [handles.smrFilePath, handles.smrFileName], 1);
if (handles.spike2file >= 0)    
    % Load Event channel 32 into memory.
    global CHAN32;
    handles.maxTime32 = cedfunction('SonChanMaxTime', handles.spike2file, 31);
    [num, CHAN32] = cedfunction('SonGetMarkData', handles.spike2file, 31, 1000000, 0, handles.maxTime32);
    handles.chand = cedfunction('SonChanDivide', handles.spike2file, 0);
    
    % Close the open spike file.
    cedfunction('SonCloseFile', handles.spike2file);
else
    disp(['Could not open SMR file: ', handles.smrFilePath, handles.smrFileName]);
    return;
end

return;



% --- Executes during object creation, after setting all properties.
function FileList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in FileList.
function FileList_Callback(hObject, eventdata, handles)
% hObject    handle to FileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns FileList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileList
