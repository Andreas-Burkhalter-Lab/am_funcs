function varargout = timestim_manuel_entry(varargin)
% TIMESTIM_MANUEL_ENTRY M-file for timestim_manuel_entry.fig
%      TIMESTIM_MANUEL_ENTRY, by itself, creates a new TIMESTIM_MANUEL_ENTRY or raises the existing
%      singleton*.
%
%      H = TIMESTIM_MANUEL_ENTRY returns the handle to a new TIMESTIM_MANUEL_ENTRY or the handle to
%      the existing singleton*.
%
%      TIMESTIM_MANUEL_ENTRY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIMESTIM_MANUEL_ENTRY.M with the given input arguments.
%
%      TIMESTIM_MANUEL_ENTRY('Property','Value',...) creates a new TIMESTIM_MANUEL_ENTRY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before timestim_manuel_entry_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to timestim_manuel_entry_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help timestim_manuel_entry

% Last Modified by GUIDE v2.5 29-Sep-2006 11:38:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @timestim_manuel_entry_OpeningFcn, ...
                   'gui_OutputFcn',  @timestim_manuel_entry_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before timestim_manuel_entry is made visible.
function timestim_manuel_entry_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to timestim_manuel_entry (see VARARGIN)

% I. Set defaults that weren't included in options

filename = varargin{1};
if size(varargin) < 2
    options = struct;
else
    options = varargin{2};
end
options = default(options,'seconds_at_a_time',60);
options = default(options,'seconds_of_overlap',20);
handles.options = options;
handles.rs.stimOnset = 0;
handles.rs.stimOffset = 0;
handles.rs.stimVolt = 0;
handles.rs.baseVolt = 0;
handles.datacursormode = 0;
handles.zoomMode = 1;
handles.dataSoFar = [];
guidata(hObject,handles)

% II. Initialize access to data:

% Open file with potential LFS support
[h,fid] = readheader(filename);

filemode = '';
is_force_use_no_lfs=ispc;
if(~is_force_use_no_lfs && strcmp(h.endian,'l'))
    if(should_use_lfs(filename))
        closelfs(fid);
    else
        fclose(fid);
    end
    [fid,msg] = openlfs(filename);
    if (fid < 0)
        error(msg);
    end
    filemode = 'lfs';
end

% Determine the stimulus channel
if ~isfield(options,'feedback')
    if(strcmp(key2value(h.wholeheader, 'hardware'), ...
            'dumb-box-0@chanel.wustl.edu'))
        options.feedback = 0;
    else
        options.feedback = 63;
    end
    warning(['warning---assuming ' num2str(options.feedback) ...
        ' is the stimulus channel in timestim().']);
end
feedback_ch_idx=find(h.channels==options.feedback);
if isempty(feedback_ch_idx)
    error(['Feedback channel ' num2str(options.feedback) ' not recorded.']);
end

% III. Load in first chunk and display

nthChnk = 1;
chnk_sz = h.scanrate * options.seconds_at_a_time;
chnk_overlap = h.scanrate * options.seconds_of_overlap;
data_start = 1;
data_end = chnk_sz + chnk_overlap;
trace = double(readthedata(data_start,data_end,h,fid,filemode,feedback_ch_idx));
nChnks = ceil(h.nscans/chnk_sz);

handles.current_data = trace;
plot(handles.current_data)
zoom on

% Put information into handles structure
handles.nthChnk = 1;
handles.chnk_sz = chnk_sz;
handles.chnk_overlap = chnk_overlap;
handles.nChnks = nChnks;
handles.xoffset = 0;
handles.h = h;
handles.fid = fid;
handles.filemode = filemode;
handles.feedback_ch_idx = feedback_ch_idx;
set(handles.progresstext,'String',['On block #' num2str(handles.nthChnk) ' of ' num2str(handles.nChnks) '!'])

% Choose default command line output for timestim_manuel_entry
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes timestim_manuel_entry wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = timestim_manuel_entry_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
handles.rs.stimOnset = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
handles.rs.stimOffset = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
handles.rs.stimVolt = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
handles.rs.baseVolt = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.datacursormode == 0
    warning('cannot accept command because data cursor mode is not active')
    return
end
info_struct =getCursorInfo(handles.dcm_obj);
pos = info_struct.Position;
if pos(1) < handles.chnk_overlap
    warning('this point is actually a part of the previous plotting region; if you need to select it, please use the "previous" button to go back')
    return
end
hold on
plotvertical(pos(1),':k')
hold off
pos(1) = pos(1) + handles.xoffset;
radiobuttonstatus = [handles.rs.stimOnset handles.rs.stimOffset handles.rs.stimVolt handles.rs.baseVolt];
handles.dataSoFar = [handles.dataSoFar; pos radiobuttonstatus];
fprintf('Data sucessfully recorded!')
guidata(hObject,handles)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nthChnk = handles.nthChnk+1;
chnk_sz = handles.chnk_sz;
chnk_overlap = handles.chnk_overlap;
data_start = (nthChnk-1)*chnk_sz+1;
data_end = min([(nthChnk*chnk_sz + chnk_overlap) handles.h.nscans]);
trace = double(readthedata(data_start,data_end,handles.h,handles.fid,handles.filemode,handles.feedback_ch_idx));

handles.current_data = trace;
hold off
plot(handles.current_data)
hold on
plotvertical(chnk_overlap,'r')
zoom on

% Put information into handles structure
handles.nthChnk = nthChnk;
handles.xoffset = data_start-1;
set(handles.progresstext,'String',['On block #' num2str(handles.nthChnk) ' of ' num2str(handles.nChnks) '!'])
set(handles.lastreminder,'Value',0)
guidata(hObject,handles)

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Shared function...
function [feedback] = readthedata(data_start,data_end,h,fid,filemode,feedback_ch_idx)

% Read the data
if strcmp(filemode,'lfs')
    tAllChData = readint16lfs(fid,h.numch,[data_start data_end], ...
        h.headersize);
else
    fseek(fid,h.numch*data_start*2+h.headersize,'bof');
    tAllChData = fread(fid,[h.numch,data_end-data_start+1],'int16');
end
feedback=tAllChData(feedback_ch_idx, :);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom on
handles.zoomMode = 1;
handles.datacursormode = 0;
guidata(hObject,handles)

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode;
handles.dcm_obj = dcm_obj;
handles.datacursormode = 1;
handles.zoomMode = 0;
guidata(hObject,handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x_y_on_off_stim_base = handles.dataSoFar;
save timestim_manuel_entry_output x_y_on_off_stim_base


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

current(1) = get(handles.radiobutton1,'Value');
current(2) = get(handles.radiobutton2,'Value');
current(3) = get(handles.radiobutton3,'Value');
current(4) = get(handles.radiobutton4,'Value');
current = (current == 0);
set(handles.radiobutton1,'Value',current(1))
set(handles.radiobutton2,'Value',current(2))
set(handles.radiobutton3,'Value',current(3))
set(handles.radiobutton4,'Value',current(4))
handles.rs.stimOffset = (handles.rs.stimOffset == 0);
handles.rs.stimOnset = (handles.rs.stimOnset == 0);
handles.rs.stimVolt = (handles.rs.stimVolt == 0);
handles.rs.baseVolt = (handles.rs.baseVolt == 0);
guidata(hObject,handles);


% --- Executes on button press in lastreminder.
function lastreminder_Callback(hObject, eventdata, handles)
% hObject    handle to lastreminder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lastreminder


