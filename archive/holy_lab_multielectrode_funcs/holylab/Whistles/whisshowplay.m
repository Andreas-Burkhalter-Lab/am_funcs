function varargout = whisshowplay(varargin)
% WHISSHOWPLAY: show and play whistles from a sng
% Syntax:
%    whisshowplay(filename,options)
%    whisshowplay(sng,header,twhis,options)
% where
%    filename is a string giving the name of the .sng file;
%    sng and header come from ReadSonogram(sngfile);
%    twhis are the whistle times (e.g., as returned by whistimes);
%    options is a structure which may have the following fields:
%      buf: gives a way to grab a bit before and after each whistle. buf
%        is a 2-vector, giving the time (in seconds) [startoffset
%        endoffset] relative to the [start end]  of each whistle.  For
%        example, if you want to see a period 5ms on either side of the
%        whistle, buf = [-0.005 0.005].
%      todisk: by  default, the play buttons produce sound from the
%        speakers.  You can instead capture the sound to disk in a file
%        called "sng.wav" (in the current directory) if todisk is true.
%      window: Or you can ignore whistles and just do the
%      window you specify as a 2 part array. 
% See also: READSONOGRAM, WHISTIMES.

% Copyright Timothy E. Holy 2004.

%WHISSHOWPLAY M-file for whisshowplay.fig
%      WHISSHOWPLAY, by itself, creates a new WHISSHOWPLAY or raises the existing
%      singleton*.
%
%      H = WHISSHOWPLAY returns the handle to a new WHISSHOWPLAY or the handle to
%      the existing singleton*.
%
%      WHISSHOWPLAY('Property','Value',...) creates a new WHISSHOWPLAY using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to whisshowplay_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      WHISSHOWPLAY('CALLBACK') and WHISSHOWPLAY('CALLBACK',hObject,...) call the
%      local function named CALLBACK in WHISSHOWPLAY.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help whisshowplay

% Last Modified by GUIDE v2.5 25-Aug-2004 16:16:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @whisshowplay_OpeningFcn, ...
    'gui_OutputFcn',  @whisshowplay_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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


% --- Executes just before whisshowplay is made visible.
function whisshowplay_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for whisshowplay
handles.output = hObject;

buf = [0 0];
options = struct;
if length(varargin)==2
    options = varargin{2};
    sng = varargin{1};
end

if (length(varargin) > 2)
    sng = varargin{1};
    header = varargin{2};
    twhis = varargin{3};
    
    if (length(varargin) > 3)
        options = varargin{4};
    end
end 
    
    
    
    
    [sng,header] = ReadSonogram(varargin{1});
    if isfield(options,'window')
        window=options.window;
        if length(window)>2
            error('just begining and end please')
        else
            tmp(1,1) =  window(1,1);
            tmp(2,1) = window(end);
            twhis = tmp;
        end
    else
        twhis = whistimes(sng,header,whistimesdefaults);
    end
    if (length(varargin) > 1)
        options = varargin{2};
    end

if isfield(options,'buf')
  buf = options.buf;
end




icons = swicons;
set(handles.leftbutton,'CData',icons.left,'Enable','off');
set(handles.rightbutton,'CData',icons.right);

handles.sng = sng;
handles.header = header;
handles.twhis = twhis;
handles.buf = buf;
handles.options = options;

set(handles.totalwhis,'String',['/' num2str(size(twhis,2))]);

wsp_display(sng,header,twhis(:,1),buf);
% Set up image properties once and for all
colormap(1-gray);
%set(gca,'YDir','reverse');
%xlabel('Time (ms)')
%ylabel('Frequency (kHz)')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes whisshowplay wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = whisshowplay_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function whisnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whisnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function shiftslowfac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shiftslowfac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Callbacks                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in playbutton. (Pitch Shift)
function playbutton_Callback(hObject, eventdata, handles)
% hObject    handle to playbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
whisnum = str2double(get(handles.whisnum,'String'));
sngsnip = wsp_display(handles.sng,handles.header,handles.twhis(:,whisnum),handles.buf);
soundsnip = sng2sound(sngsnip);
shiftfac = str2double(get(handles.shiftslowfac,'String'));
freq = [handles.header.freqMin handles.header.freqMax]/handles.header.nfreq;
shiftsound = phasevocoder(soundsnip,2*(handles.header.nfreq-1),shiftfac,handles.header.threshold,freq);
mass = max(abs(shiftsound));
if (mass > 0)
  shiftsound = shiftsound/(1.001*mass);
end
if (isfield(handles.options,'todisk') & handles.options.todisk)
  wavwrite(shiftsound,handles.header.scanrate/shiftfac,'sng.wav');
else
  sound(shiftsound,handles.header.scanrate/shiftfac);
end


% --- Executes on button press in slowdownbutton.
function slowdownbutton_Callback(hObject, eventdata, handles)
% hObject    handle to slowdownbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
whisnum = str2double(get(handles.whisnum,'String'));
sngsnip = wsp_display(handles.sng,handles.header,handles.twhis(:,whisnum),handles.buf);
soundsnip = sng2sound(sngsnip);
slowfac = str2double(get(handles.shiftslowfac,'String'));
mass = max(abs(soundsnip));
if (mass > 0)
  soundsnip = soundsnip/(1.001*mass);
end
if (isfield(handles.options,'todisk') & handles.options.todisk)
  wavwrite(soundsnip,handles.header.scanrate/slowfac,'sng.wav');
else
  sound(soundsnip,handles.header.scanrate/slowfac);
end


function shiftslowfac_Callback(hObject, eventdata, handles)
% hObject    handle to shiftslowfac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shiftslowfac as text
%        str2double(get(hObject,'String')) returns contents of shiftslowfac as a double

% No action required





% --- Executes on button press in leftbutton.
function leftbutton_Callback(hObject, eventdata, handles)
% hObject    handle to leftbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
whisnum = str2double(get(handles.whisnum,'String'));
if (whisnum > 1)
  whisnum = whisnum-1;
  set(handles.whisnum,'String',int2str(whisnum));
end
wsp_display(handles.sng,handles.header,handles.twhis(:,whisnum),handles.buf);
wsp_enable(handles,whisnum);

% --- Executes on button press in rightbutton.
function rightbutton_Callback(hObject, eventdata, handles)
% hObject    handle to rightbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
whisnum = str2double(get(handles.whisnum,'String'));
if (whisnum < size(handles.twhis,2))
  whisnum = whisnum+1;
  set(handles.whisnum,'String',int2str(whisnum));
end
wsp_display(handles.sng,handles.header,handles.twhis(:,whisnum),handles.buf);
wsp_enable(handles,whisnum);

function whisnum_Callback(hObject, eventdata, handles)
% hObject    handle to whisnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whisnum as text
%        str2double(get(hObject,'String')) returns contents of whisnum as a double
whisnum = str2double(get(handles.whisnum,'String'));
if (whisnum < 1)
  whisnum = 1;
elseif (whisnum > size(handles.twhis,2))
  whisnum = size(handles.twhis,2);
end
set(handles.whisnum,'String',int2str(whisnum));
wsp_display(handles.sng,handles.header,handles.twhis(:,whisnum),handles.buf);
wsp_enable(handles,whisnum);


function ws = wsp_display(sng,header,t,buf)
T = header.nscans/header.scanrate;
t1 = t + buf';
f = [0 header.scanrate/2000];
tc = round(t1 * (header.columnTotal - 1)/T + 1);
ws = sng(:,tc(1):tc(2));
if (nargout == 0)
  imagesc(t1,f,abs(ws));
  line((t*[1 1])',[f' f'],'Color','r');
  set(gca,'YDir','normal');
  xlabel('Time (s)')
  ylabel('Frequency (kHz)')
end

function wsp_enable(handles,whisnum)
if (whisnum > 1)
  set(handles.leftbutton,'Enable','on')
else
  set(handles.leftbutton,'Enable','off')
end
if (whisnum < size(handles.twhis,2))
  set(handles.rightbutton,'Enable','on')
else
  set(handles.rightbutton,'Enable','off')
end


% A helper function to draw the paging icons
% Copied from sliderwindow
function icons = swicons
  iconsize = [17 17];
  bkgrnd = get(0,'DefaultUicontrolBackgroundColor');
  blank = zeros(iconsize) + bkgrnd(1);
  
  right = blank;
  index = mask(iconsize,1,1);
  right(index) = 0;
  right = repmat(right,[1 1 3]);

  rightsm = blank;
  index = mask(iconsize,0.5,1);
  rightsm(index) = 0;
  rightsm = repmat(rightsm,[1 1 3]);

  left = blank;
  index = mask(iconsize,1,-1);
  left(index) = 0;
  left = repmat(left,[1 1 3]);
  
  leftsm = blank;
  index = mask(iconsize,0.5,-1);
  leftsm(index) = 0;
  leftsm = repmat(leftsm,[1 1 3]);
  
  icons.right = right;
  icons.rightsm = rightsm;
  icons.left = left;
  icons.leftsm = leftsm;
  
% This function takes parameters from swicons
% to do the actual icon "drawing"
function index = mask(iconsize,width,direction)
  [row,col] = find(ones(iconsize));
  % Make them go from 0 to 1
  row = (row-1)/(iconsize(1)-1);
  col = (col-1)/(iconsize(2)-1);
  if (direction > 0)
    index = find(2*abs(row-0.5) <= (width-col)/width);
  else
    index = find(2*abs(row-0.5) + (1 - width)/width <= col/width);
  end
  return;


