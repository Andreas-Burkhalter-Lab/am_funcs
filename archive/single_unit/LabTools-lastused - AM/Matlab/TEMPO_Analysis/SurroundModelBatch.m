function varargout = SurroundModelBatch(varargin)
% SURROUNDMODELBATCH Application M-file for SurroundModelBatch.fig
%    FIG = SURROUNDMODELBATCH launch SurroundModelBatch GUI.
%    SURROUNDMODELBATCH('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 16-Oct-2001 13:54:42

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
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
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
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
    PATH = get(handles.PathEdit, 'String');
    t = sprintf('%s*.m', PATH);
    [fname, pname] = uigetfile(t,'Select batch file');	
    if (fname ~= 0)	
        set(handles.PathEdit, 'String', pname);
        set(handles.FileEdit, 'String', fname);
    end



% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton2.
disp('WA-HOO!!!')

   batchPATH = get(handles.PathEdit, 'String');
   batchFILE = get(handles.FileEdit, 'String');
   filename = [batchPATH batchFILE];
   fid = fopen(filename);

   line = fgetl(fid);
   %1) while we are not at the EOF do:   
   while(line ~= -1)
      %skip any comment lines
      while(line(1) == '%')
         line = fgetl(fid);
      end
      
      %  2) we want to read in the next 3 lines in the file
      
      %first remove any comment text at the end of a line (following a %)
      comment_start = find(line == '%');
      if ~isempty(comment_start)
         line = line(1:(comment_start(1)-1));
      end
            
      %get an index of spaces so we know where data starts and stops
      spaces = isspace(line);
      space_index = find(spaces);
        
      %get path / file
      %     3) the first line is the disparity file
      PATH = line(1:space_index(1) - 1);
      DispFILE = line(space_index(1) + 1:space_index(2) - 1)
            
      DispBegTrial = line(space_index(2) + 1:space_index(3) - 1);
      DispBegTrial = str2num(DispBegTrial);
      
      DispEndTrial = line(space_index(3) + 1:space_index(4) - 1);
      DispEndTrial = str2num(DispEndTrial);
      
      DispStartCode = line(space_index(4) + 1:space_index(5) - 1);
      DispStartCode = str2num(DispStartCode);
      
      DispStartOffset = line(space_index(5) + 1:space_index(6) - 1);
      DispStartOffset = str2num(DispStartOffset);
      
      DispStopCode = line(space_index(6) + 1:space_index(7) - 1);
      DispStopCode = str2num(DispStopCode);
      
      DispStopOffset = line(space_index(7) + 1:space_index(8) - 1);
      DispStopOffset = str2num(DispStopOffset);
      
      DispSpikeChan = line(space_index(8) + 1:space_index(9) - 1);
      DispSpikeChan = str2num(DispSpikeChan);
      
      if length(line) == space_index(9) + 1
         good_flag = line(space_index(9)+1:length(line));
         good_flag = str2num(good_flag);
      else
         good_flag = line(space_index(9) + 1:space_index(10) - 1);
         good_flag = str2num(good_flag);
      end
            
      %first load hdisp data
      %     4) load the disparity file
      [HDisp_return_val, HDisp_good_data, HDisp_bad_data] = LoadTEMPOData(PATH,DispFILE);

            
      if (good_flag)
         Disp_select_data = HDisp_good_data;
      else
         Disp_select_data = HDisp_bad_data;
      end
      
      if DispBegTrial == -1
         DispBegTrial = 1;		
      end
            
      if DispEndTrial == -1
         DispEndTrial = size(Disp_select_data.event_data,3);		
      end
      
      %     5) the second line is the gradient file
      %skip any comment lines that might exist before the next data file
      line = fgetl(fid);
      while(line(1) == '%')
         line = fgetl(fid);
      end
      
      %first remove any comment text at the end of a line (following a %)
      comment_start = find(line == '%');
      if ~isempty(comment_start)
         line = line(1:(comment_start(1)-1));
      end
            
      %get an index of spaces so we know where data starts and stops
      spaces = isspace(line);
      space_index = find(spaces);
        
      %get path / file
      %     3) the first line is the disparity file
      PATH = line(1:space_index(1) - 1);
      GradFILE = line(space_index(1) + 1:space_index(2) - 1)
      
      GradBegTrial = line(space_index(2) + 1:space_index(3) - 1);
      GradBegTrial = str2num(GradBegTrial);
      
      GradEndTrial = line(space_index(3) + 1:space_index(4) - 1);
      GradEndTrial = str2num(GradEndTrial);
      
      GradStartCode = line(space_index(4) + 1:space_index(5) - 1);
      GradStartCode = str2num(GradStartCode);
      
      GradStartOffset = line(space_index(5) + 1:space_index(6) - 1);
      GradStartOffset = str2num(GradStartOffset);
      
      GradStopCode = line(space_index(6) + 1:space_index(7) - 1);
      GradStopCode = str2num(GradStopCode);
      
      GradStopOffset = line(space_index(7) + 1:space_index(8) - 1);
      GradStopOffset = str2num(GradStopOffset);
    
      GradSpikeChan = line(space_index(8) + 1:space_index(9) - 1);
      GradSpikeChan = str2num(GradSpikeChan);

      %depending on how the Good Trials/Bad Trials radio buttons are set, analyze the appropriate set of trials
      if length(line) == space_index(9) + 1
         good_flag = line(space_index(9)+1:length(line));
         good_flag = str2num(good_flag);
      else
         good_flag = line(space_index(9) + 1:space_index(10) - 1);
         good_flag = str2num(good_flag);
      end
            
      %     6) load the gradient file
      [Grad_return_val, Grad_good_data, Grad_bad_data] = LoadTEMPOData(PATH,GradFILE);
     
      
      if (good_flag)
          Grad_select_data = Grad_good_data;
      else
          Grad_select_data = Grad_bad_data;
      end
      
      if GradBegTrial == -1
         GradBegTrial = 1;		
      end
            
      if GradEndTrial == -1
         GradEndTrial = size(Grad_select_data.event_data,3);		
      end
      
      %     7) the third line is the surround mapping file
      %     8) load the sur. map file
      
      line = fgetl(fid);
      %skip any comment lines that might exist before the next data file
      while(line(1) == '%')
         line = fgetl(fid);
      end
      
      %first remove any comment text at the end of a line (following a %)
      comment_start = find(line == '%');
      if ~isempty(comment_start)
         line = line(1:(comment_start(1)-1));
      end
            
      %get an index of spaces so we know where data starts and stops
      spaces = isspace(line);
      space_index = find(spaces);
        
      %get path / file
      %     3) the first line is the disparity file
      PATH = line(1:space_index(1) - 1);
      SurrMapFILE = line(space_index(1) + 1:space_index(2) - 1)
      
      SurrMapBegTrial = line(space_index(2) + 1:space_index(3) - 1);
      SurrMapBegTrial = str2num(SurrMapBegTrial);
      
      SurrMapEndTrial = line(space_index(3) + 1:space_index(4) - 1);
      SurrMapEndTrial = str2num(SurrMapEndTrial);
      
      SurrMapStartCode = line(space_index(4) + 1:space_index(5) - 1);
      SurrMapStartCode = str2num(SurrMapStartCode);
      
      SurrMapStartOffset = line(space_index(5) + 1:space_index(6) - 1);
      SurrMapStartOffset = str2num(SurrMapStartOffset);
      
      SurrMapStopCode = line(space_index(6) + 1:space_index(7) - 1);
      SurrMapStopCode = str2num(SurrMapStopCode);
      
      SurrMapStopOffset = line(space_index(7) + 1:space_index(8) - 1);
      SurrMapStopOffset = str2num(SurrMapStopOffset);
    
      SurrMapSpikeChan = line(space_index(8) + 1:space_index(9) - 1);
      SurrMapSpikeChan = str2num(SurrMapSpikeChan);

      %depending on how the Good Trials/Bad Trials radio buttons are set, analyze the appropriate set of trials
      if length(line) == space_index(9) + 1
         good_flag = line(space_index(9)+1:length(line));
         good_flag = str2num(good_flag);
      else
         good_flag = line(space_index(9) + 1:space_index(10) - 1);
         good_flag = str2num(good_flag);
      end
            
      %     6) load the gradient file
      [SurrMap_return_val, SurrMap_good_data, SurrMap_bad_data] = LoadTEMPOData(PATH,SurrMapFILE);
      
      if (good_flag)
          SurrMap_select_data = SurrMap_good_data;
      else
          SurrMap_select_data = SurrMap_bad_data;
      end
      
      if SurrMapBegTrial == -1
         SurrMapBegTrial = 1;		
      end
            
      if SurrMapEndTrial == -1
         SurrMapEndTrial = size(SurrMap_select_data.event_data,3);		
      end
      
      %now that we have loaded all the data, we need to analyze the heck out of it!!  ANALYZE GO!!
      UseSyncPulses = 0;            
      
      %calculate the spike rates for the each data file
      %ok, now what?  we need to process the disparity data and obtain the parameters of the tuning curve
      HDisparity.perr = [];
      HDisparity.px = [];
      HDisparity.py = [];
      HDisparity.uncorr_resp = [];
      HDisparity.plot_x = [];
      HDisparity.plot_y = [];

      [HDisparity.perr, HDisparity.px, HDisparity.py, HDisparity.uncorr_resp, HDisparity.plot_x, HDisparity.plot_y] = HDispTuningParams(Disp_select_data, DispSpikeChan, DispStartCode, DispStopCode, DispBegTrial, DispEndTrial, DispStartOffset, DispStopOffset, PATH, DispFILE, UseSyncPulses);
   
      %for surround mapping i need to obtain the following info:
      %surround size
      %center size

      %and for each surround patch:
         %perr
         %px
         %py
         %uncorr_resp
         %plot_x
         %plot_y
         
      %now we need the parameters of the tuning curves of each of the surround patches
      SurrMap.ctronly = [];
      SurrMap.ctrdisp = [];
      SurrMap.surrsize = [];
      SurrMap.ctrsize = [];
      SurrMap.angles = [];
      SurrMap.px = [];
      SurrMap.py = [];
      SurrMap.plot_x = [];
      SurrMap.plot_y = [];
      [SurrMap.ctronly, SurrMap.ctrdisp, SurrMap.surrsize, SurrMap.ctrsize, SurrMap.angles, SurrMap.px, SurrMap.py, SurrMap.plot_x, SurrMap.plot_y] = SurrMapTuningParams(SurrMap_select_data, SurrMapSpikeChan, SurrMapStartCode, SurrMapStopCode, SurrMapBegTrial, SurrMapEndTrial, SurrMapStartOffset, SurrMapStopOffset, PATH, SurrMapFILE, UseSyncPulses);
      
      %now we can simulate the response of the neuron to h. disparity gradients at different orientations
      SurroundModel2(Grad_select_data, GradSpikeChan, GradStartCode, GradStopCode, GradBegTrial, GradEndTrial, GradStartOffset, GradStopOffset, PATH, GradFILE, HDisparity, SurrMap, UseSyncPulses);
      
      print_flag = get(handles.PrintHandle, 'Value');
      if print_flag
         printhandle = figure;
         close(printhandle);
         for printindex = 1:1:printhandle - 1
            print(printindex);
         end
      end
         
      close_flag = get(handles.CloseHandle, 'Value');
      if close_flag
         closehandle = figure;
         for closeindex = 1:1:closehandle
            close(closeindex);
         end
      end
      
      
      %get next line
      line = fgetl(fid);

   end %end while for EOF
            
            
% end function pushbutton2_Callback----------------------------------------------------------
