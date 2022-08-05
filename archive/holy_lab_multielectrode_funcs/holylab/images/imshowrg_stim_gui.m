function varargout = imshowrg_stim_gui(varargin)
% Syntax Example:
% imshowrg_stim_gui('vno_gcamp2_2009_12_05.imagine')
%
% In case imagine crashed in the middle of acquisition, you have to go and
% change the .imagine file's stack number field to represent the last
% working stack. Then this gui will work.
%
% works for upto 15 unique stimuli (with one flush). upto 5 repeats
% plots activity as red/green according to the function colorize_dfof
% takes the average of 4 frames before stimulus and average of 4 frames
% after stimulus and calculated dfof image

% IMSHOWRG_STIM_GUI M-file for imshowrg_stim_gui.fig
%      IMSHOWRG_STIM_GUI, by itself, creates a new IMSHOWRG_STIM_GUI or raises the existing
%      singleton*.
%
%      H = IMSHOWRG_STIM_GUI returns the handle to a new IMSHOWRG_STIM_GUI or the handle to
%      the existing singleton*.
%
%      IMSHOWRG_STIM_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMSHOWRG_STIM_GUI.M with the given input arguments.
%
%      IMSHOWRG_STIM_GUI('Property','Value',...) creates a new IMSHOWRG_STIM_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imshowrg_stim_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imshowrg_stim_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% Last Modified by GUIDE v2.5 13-Jan-2010 21:06:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imshowrg_stim_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @imshowrg_stim_gui_OutputFcn, ...
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


% --- Executes just before imshowrg_stim_gui is made visible.
function imshowrg_stim_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imshowrg_stim_gui (see VARARGIN)

% Choose default command line output for imshowrg_stim_gui
handles.output = hObject;

filename = varargin{1};

smm = stackmm(filename);

sz = smm.size;
header = smm.header;
stim = header.stim_lookup;
stim_labels = header.stim_labels;

flush = mode(stim); % flush is the most common valve
stim_id = unique(stim); % all valves
stim_id = setxor(flush, stim_id); % all stimuli...now with flush removed

stim_start = nan(16, 5); % 16 valves max, 5 repeats max

for stim_indx = 1:length(stim_id)
    count = 1;
    for indx = 2:sz(4)
        if stim(indx-1)==flush && stim(indx)==stim_id(stim_indx)
            stim_start(stim_id(stim_indx), count) = indx;
            count = count+1;
        end
    end
end

setappdata(handles.figStimGUI,'smm',smm);
setappdata(handles.figStimGUI,'stim_labels', stim_labels);
setappdata(handles.figStimGUI,'stim_id',stim_id);
setappdata(handles.figStimGUI,'stim_start',stim_start);

clim(1) = str2num(get(handles.clim_min,'String'));
clim(2) = str2num(get(handles.clim_max,'String'));

dfof(1) = str2num(get(handles.dfof_neg,'String'));
dfof(2) = str2num(get(handles.dfof_pos,'String'));

frame_num = str2num(get(handles.frame_number,'String'));
stim = get(handles.current_stimulus, 'String');
stim = stim_labels{stim_id(1)}; % in case stim_label1 = flush
stim_indx = stim_id(1); % in case there are multiple stimuli with same name

setappdata(handles.figStimGUI,'clim',clim);
setappdata(handles.figStimGUI,'dfof',dfof);
setappdata(handles.figStimGUI,'frame_num',frame_num);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);

set_labels_on_buttons(handles);

plot_frames(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imshowrg_stim_gui wait for user response (see UIRESUME)
% uiwait(handles.figStimGUI);


% --- Outputs from this function are returned to the command line.
function varargout = imshowrg_stim_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in show_stim1.
function show_stim1_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{1};
stim_indx = 1;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);

% --- Executes on button press in show_stim2.
function show_stim2_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{2};
stim_indx = 2;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);


% --- Executes on button press in show_stim3.
function show_stim3_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
%setappdata(handles.figStimGUI,'stim','Ringers');
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{3};
stim_indx = 3;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);

% --- Executes on button press in show_stim4.
function show_stim4_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{4};
stim_indx = 4;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);



% --- Executes on button press in show_stim5.
function show_stim5_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{5};
stim_indx = 5;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);


% --- Executes on button press in show_stim6.
function show_stim6_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{6};
stim_indx = 6;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);


% --- Executes on button press in show_stim7.
function show_stim7_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{7};
stim_indx = 7;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);



% --- Executes on button press in show_stim8.
function show_stim8_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{8};
stim_indx = 8;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);



% --- Executes on button press in show_stim9.
function show_stim9_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{9};
stim_indx = 9;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);



% --- Executes on button press in show_stim10.
function show_stim10_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{10};
stim_indx = 10;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);



% --- Executes on button press in show_stim11.
function show_stim11_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{11};
stim_indx = 11;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);



% --- Executes on button press in show_stim12.
function show_stim12_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{12};
stim_indx = 12;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);



% --- Executes on button press in show_stim13.
function show_stim13_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{13};
stim_indx = 13;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);



% --- Executes on button press in show_stim14.
function show_stim14_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{14};
stim_indx = 14;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);

% --- Executes on button press in show_stim15.
function show_stim15_Callback(hObject, eventdata, handles)
% hObject    handle to show_stim15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim = stim_labels{15};
stim_indx = 15;
set(handles.current_stimulus, 'String', stim);
setappdata(handles.figStimGUI,'stim',stim);
setappdata(handles.figStimGUI,'stim_indx',stim_indx);
plot_frames(handles);


function clim_min_Callback(hObject, eventdata, handles)
% hObject    handle to clim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clim_min as text
%        str2double(get(hObject,'String')) returns contents of clim_min as a double
clim(1) = str2num(get(handles.clim_min,'String'));
clim(2) = str2num(get(handles.clim_max,'String'));
setappdata(handles.figStimGUI,'clim',clim);
plot_frames(handles);

% --- Executes during object creation, after setting all properties.
function clim_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clim_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function clim_max_Callback(hObject, eventdata, handles)
% hObject    handle to clim_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clim_max as text
%        str2double(get(hObject,'String')) returns contents of clim_max as a double
clim(1) = str2num(get(handles.clim_min,'String'));
clim(2) = str2num(get(handles.clim_max,'String'));
setappdata(handles.figStimGUI,'clim',clim);
plot_frames(handles);

% --- Executes during object creation, after setting all properties.
function clim_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clim_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dfof_neg_Callback(hObject, eventdata, handles)
% hObject    handle to dfof_neg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dfof_neg as text
%        str2double(get(hObject,'String')) returns contents of dfof_neg as a double
dfof(1) = str2num(get(handles.dfof_neg,'String'));
dfof(2) = str2num(get(handles.dfof_pos,'String'));
setappdata(handles.figStimGUI,'dfof',dfof);
plot_frames(handles);


% --- Executes during object creation, after setting all properties.
function dfof_neg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dfof_neg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dfof_pos_Callback(hObject, eventdata, handles)
% hObject    handle to dfof_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dfof_pos as text
%        str2double(get(hObject,'String')) returns contents of dfof_pos as a double
dfof(1) = str2num(get(handles.dfof_neg,'String'));
dfof(2) = str2num(get(handles.dfof_pos,'String'));
setappdata(handles.figStimGUI,'dfof',dfof);
plot_frames(handles);


% --- Executes during object creation, after setting all properties.
function dfof_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dfof_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_number_Callback(hObject, eventdata, handles)
% hObject    handle to frame_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_number as text
%        str2double(get(hObject,'String')) returns contents of frame_number as a double
frame_num = str2num(get(handles.frame_number,'String'));
setappdata(handles.figStimGUI,'frame_num',frame_num);
plot_frames(handles);


% --- Executes during object creation, after setting all properties.
function frame_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in prev_frame.
function prev_frame_Callback(hObject, eventdata, handles)
% hObject    handle to prev_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frame_num = str2num(get(handles.frame_number,'String'));
frame_num = frame_num - 1;
if frame_num > 0
    set(handles.frame_number,'String',num2str(frame_num));
    setappdata(handles.figStimGUI,'frame_num',frame_num);
end
plot_frames(handles);
    

% --- Executes on button press in next_frame.
function next_frame_Callback(hObject, eventdata, handles)
% hObject    handle to next_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

frame_num = str2num(get(handles.frame_number,'String'));
smm = getappdata(handles.figStimGUI,'smm');
sz = smm.size;

frame_num = frame_num + 1;
if frame_num <= sz(3)
    set(handles.frame_number,'String',num2str(frame_num));
    setappdata(handles.figStimGUI,'frame_num',frame_num);
end
plot_frames(handles);

function plot_frames(handles)

smm = getappdata(handles.figStimGUI,'smm');
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim_id = getappdata(handles.figStimGUI,'stim_id');
stim_start = getappdata(handles.figStimGUI,'stim_start');
clim = getappdata(handles.figStimGUI,'clim');
dfof = getappdata(handles.figStimGUI,'dfof');
frame_num = getappdata(handles.figStimGUI,'frame_num');
stim = getappdata(handles.figStimGUI,'stim');
sz = smm.size;

%current_stim_loc = findainb(stim, stim_labels);
current_stim_loc = getappdata(handles.figStimGUI,'stim_indx');

for axes_indx = 1:5

    indx = stim_start(current_stim_loc, axes_indx);

    if ~isnan(indx) % making sure the stimulus repeat exists

        im1 = mean(smm(:,:,frame_num, indx-5:indx-1), 4); %bg
        im2 = mean(smm(:,:,frame_num, indx:indx+4), 4); %stim

        im(:,:,1) = im1;
        im(:,:,2) = im2;
        
        frames{1} = im;
        options.clim_dfof = dfof;
        %options.sigma = 2;
        options.clim_raw = clim;

        imrgb_tmp = colorize_dfof(frames, options);
        imrgb = imrotate(imrgb_tmp{1},90);
        
    else

        imrgb = zeros([sz(1) sz(2)]);
    end

    if axes_indx == 1
        axes(handles.axes1);
    elseif axes_indx == 2
        axes(handles.axes2);
    elseif axes_indx == 3
        axes(handles.axes3);
    elseif axes_indx == 4
        axes(handles.axes4);
    elseif axes_indx == 5
        axes(handles.axes5);
    end




    imagesc(imrgb);
    axis image;
    axis off;



end

function set_labels_on_buttons(handles)
stim_labels = getappdata(handles.figStimGUI,'stim_labels');
stim_id = getappdata(handles.figStimGUI, 'stim_id');

if ~isempty(intersect(stim_id, 1))
    set(handles.show_stim1,'String',stim_labels{1})
else
    delete(handles.show_stim1);
end

if ~isempty(intersect(stim_id, 2))
    set(handles.show_stim2,'String',stim_labels{2})
else
    delete(handles.show_stim2);
end

if ~isempty(intersect(stim_id, 3))
    set(handles.show_stim3,'String',stim_labels{3})
else
    delete(handles.show_stim3);
end

if ~isempty(intersect(stim_id, 4))
    set(handles.show_stim4,'String',stim_labels{4})
else
    delete(handles.show_stim4);
end

if ~isempty(intersect(stim_id, 5))
    set(handles.show_stim5,'String',stim_labels{5})
else
    delete(handles.show_stim5);
end

if ~isempty(intersect(stim_id, 6))
    set(handles.show_stim6,'String',stim_labels{6})
else
    delete(handles.show_stim6);
end

if ~isempty(intersect(stim_id, 7))
    set(handles.show_stim7,'String',stim_labels{7})
else
    delete(handles.show_stim7);
end

if ~isempty(intersect(stim_id, 8))
    set(handles.show_stim8,'String',stim_labels{8})
else
    delete(handles.show_stim8);
end

if ~isempty(intersect(stim_id, 9))
    set(handles.show_stim9,'String',stim_labels{9})
else
    delete(handles.show_stim9);
end

if ~isempty(intersect(stim_id, 10))
    set(handles.show_stim10,'String',stim_labels{10})
else
    delete(handles.show_stim10);
end

if ~isempty(intersect(stim_id, 11))
    set(handles.show_stim11,'String',stim_labels{11})
else
    delete(handles.show_stim11);
end

if ~isempty(intersect(stim_id, 12))
    set(handles.show_stim12,'String',stim_labels{12})
else
    delete(handles.show_stim12);
end

if ~isempty(intersect(stim_id, 13))
    set(handles.show_stim13,'String',stim_labels{13})
else
    delete(handles.show_stim13);
end

if ~isempty(intersect(stim_id, 14))
    set(handles.show_stim14,'String',stim_labels{14})
else
    delete(handles.show_stim14);
end

if ~isempty(intersect(stim_id, 15))
    set(handles.show_stim15,'String',stim_labels{15})
else
    delete(handles.show_stim15);
end





