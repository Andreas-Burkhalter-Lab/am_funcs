function varargout = rec_whis2(varargin)
% REC_WHIS2 M-file for rec_whis2.fig
%      REC_WHIS2, by itself, creates a new REC_WHIS2 or raises the existing
%      singleton*.
%
%      H = REC_WHIS2 returns the handle to a new REC_WHIS2 or the handle to
%      the existing singleton*.
%
%      REC_WHIS2('Property','Value',...) creates a new REC_WHIS2 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to rec_whis2_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      REC_WHIS2('CALLBACK') and REC_WHIS2('CALLBACK',hObject,...) call the
%      local function named CALLBACK in REC_WHIS2.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
   
% Edit the above text to modify the response to help rec_whis2
   
% Last Modified by GUIDE v2.5 25-Jul-2003 11:47:50
   
% Begin initialization code - DO NOT EDIT
   gui_Singleton = 1;
   gui_State = struct('gui_Name',       mfilename, ...
                      'gui_Singleton',  gui_Singleton, ...
                      'gui_OpeningFcn', @rec_whis2_OpeningFcn, ...
                      'gui_OutputFcn',  @rec_whis2_OutputFcn, ...
                      'gui_LayoutFcn',  [], ...
                      'gui_Callback',   []);
   if nargin & isstr(varargin{1})
      gui_State.gui_Callback = str2func(varargin{1});
   end
   
   if nargout
      [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
   else
      gui_mainfcn(gui_State, varargin{:});
   end
% End initialization code - DO NOT EDIT


% --- Executes just before rec_whis2 is made visible.
function rec_whis2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
   
% Choose default command line output for rec_whis2
   handles.output = hObject;
   
   % Update handles structure
   guidata(hObject, handles);
   
   % UIWAIT makes rec_whis2 wait for user response (see UIRESUME)
   % uiwait(handles.figure1);
   
   editSongScanRate_Callback(hObject, [], handles);
   
   
% --- Outputs from this function are returned to the command line.
function varargout = rec_whis2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
% Get default command line output from handles structure
   varargout{1} = handles.output;
   
   
% --- Executes on button press in btnTest.
function btnTest_Callback(hObject, eventdata, handles)
% hObject    handle to btnTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   msgbox('hello, new guide');
   
   
% --- Executes on button press in btnCloseMain.
function btnCloseMain_Callback(hObject, eventdata, handles)
% hObject    handle to btnCloseMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   close(handles.RecWhisParamsFig);
   
   
% ----------------------------------------------------------------
% shared by several callbacks:
function [oScanrate, oRatio, o_nFreq,o_nAvg,o_nDispColPerTransfer,o_nTransfers,oRealWholeDuration, oSongChan, oRangeIndex, ...
         o_nScansPerTransfer] = RWCBGetData(handles)
   songScanrate  = str2num(get(handles.editSongScanRate,'String'));
   sensorScanrate= str2num(get(handles.editSensorScanRate,'String')); 
   oScanrate=sensorScanrate; % sensor's scanrate is same as channel list's scanrate
   oRatio=round(songScanrate/sensorScanrate); % the channel list will be:
                                              % oRatio-1 # of song_channel and 1 sensor_channel
   songScanrate=oRatio*sensorScanrate; % the ratio must be an integer
   set(handles.editSongScanRate, 'String', num2str(songScanrate));
   
   o_nFreq = str2num(get(handles.editNumOfFreq,'String'));
   o_nDispColPerTransfer = str2num(get(handles.nDispColPerTransfer,'String')); % # of display columns per transfer
   o_nAvg = str2num(get(handles.Navg,'String')); % # of FFT columns used by each display column
   tnFFTColPerTransfer=o_nDispColPerTransfer*o_nAvg; % # of FFT columns per transfer
   tnSamplePerTransfer=tnFFTColPerTransfer*(2*o_nFreq); % # of samples per transfer
   o_nScansPerTransfer=ceil(tnSamplePerTransfer/oRatio); % # of scans per transfer
   
   tTransferInterval=o_nScansPerTransfer/oScanrate; % transfer interval in seconds
   
   acqtime = str2num(get(handles.AcqTime,'String'));
   o_nTransfers = round(acqtime/tTransferInterval); % # of transfers, i.e. # of transfer intervals
   oRealWholeDuration = o_nTransfers*tTransferInterval;

   set(handles.Info2,'String',sprintf('Transfer interval: %1.2fs\nActual acquisition time: %4.2fs', tTransferInterval, oRealWholeDuration));
   
   oSongChan = get(handles.comboxSongChannel,'Value') - 1;
   oRangeIndex= get(handles.comboxRangeIndex,'Value') - 1;

   return
   
   
% --- Executes on button press in OKButton.
% @notes: previously it is "RecWhisCB OKButton"
function OKButton_Callback(hObject, eventdata, handles)
% hObject    handle to OKButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   hfig=handles.RecWhisParamsFig;
   [p.scanrate, p.ratio, p.nFreq,p.nAvg,p.nDispColPerTransfer,p.nTransfers,p.realWholeDuration,p.songchan,p.rangeIndex, ...
   p.nScansPerTransfer] = RWCBGetData(handles);
   if(get(handles.checkboxSaveSongRawData,'Value')) % if save song's raw data
      p.proxDetect = get(handles.checkboxProxDetect,'Value');
      if (p.proxDetect)
         p.sensorchan = get(handles.comboxSensorChan,'Value') ;
      else
         p.sensorchan = [ ];
      end    
      if (p.songchan == p.sensorchan)
         errordlg('Channels cannot be the same','Input Error');
         return
      end
      temphdr = get(handles.UsrHdr,'String')';
      ct = cellstr(temphdr');                % Pad header with returns at end of lines
      ret = double(sprintf('\n'));
      ctcat = strcat(ct,num2cell(char(ret*ones(length(ct),1))));
      p.usrhdr = cat(2,ctcat{:})';        % Turn into a single column string
      [fname,pathname] = uiputfile('*.bin','Save data to file:');
      if (fname == 0)
         disp('Operation cancelled');
         return
      end
      p.filename = [pathname,fname];
      p.savefile = 1;
   else % else, not save
      p.usrhdr = [];
      p.filename = [];
      p.savefile = 0;
      p.sensorchan = [ ];
      p.proxDetect = 0;
   end
   assignin('base','recparams',p);
   % tthRW=DoRW(p); % get handle to be used by guide(tthRW)
   tthRW=do_rw2(p);
   ttInt=9; % just provide a breakpoint for debug
   
   
   
% ---------------------------------------------------------------
% @notes: previously it is called "RecWhisCB UpdateInfo";
%    shared by editSongScanRate, AcqTime, editNumOfFreq, Navg, nDispColPerTransfer, checkboxSaveSongRawData, checkboxProxDetect:
function editSongScanRate_Callback(hObject, eventdata, handles)
% hObject    handle to editSongScanRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
% Hints: get(hObject,'String') returns contents of editSongScanRate as text
%        str2double(get(hObject,'String')) returns contents of editSongScanRate as a double
   hfig=handles.RecWhisParamsFig;
   [scanrate, ratio, nFreq,nAvg,nDispColPerTransfer,nTransfers,realWholeDuration,songchan,rangeIndex, nScansPerTransfer] = RWCBGetData(handles);
   
   if(get(handles.checkboxSaveSongRawData,'Value'))
      set(handles.UsrHdr,'Enable','on');
      set(handles.checkboxProxDetect,'Value',1, 'enable', 'on');
      set(handles.comboxSensorChan,'Enable','on');
      set(handles.StaticText8,'Enable','on');
   else 
      set(handles.UsrHdr,'Enable','off');
      set(handles.checkboxProxDetect,'Value',0, 'Enable','off');
      set(handles.comboxSensorChan,'Enable','off');
      set(handles.StaticText8,'Enable','off');
   end
   %------------------------------------------------------------------------------
   


% --- Executes during object creation, after setting all properties.
function editSensorScanRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSensorScanRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
return;
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes during object creation, after setting all properties.
function editSongScanRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSongScanRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
return;
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


