function varargout = SurroundModel_GUI(varargin)
% SURROUNDMODEL_GUI Application M-file for SurroundModel_GUI.fig
%    FIG = SURROUNDMODEL_GUI launch SurroundModel_GUI GUI.
%    SURROUNDMODEL_GUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 12-Mar-2001 14:44:56

global HDisp_good_data;
global SurrMap_good_data
global Grad_good_data
global HDisp_bad_data;
global SurrMap_bad_data;
global Grad_bad_data;
global listtext;

load SurroundModel_GUI
TEMPO_Defs;
Path_Defs;

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
   set(fig, 'Name', 'Surround Model GUI');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
   
   %setup popupmenu items
   set(handles.DispStartTimeMenu, 'String', mat1);
   set(handles.DispStartTimeMenu, 'Value', 4);
   set(handles.DispStopTimeMenu, 'String', mat2);
   set(handles.DispStopTimeMenu, 'Value', 5);
   
   set(handles.SurrMapStartTimeMenu, 'String', mat1);
   set(handles.SurrMapStartTimeMenu, 'Value', 4);
   set(handles.SurrMapStopTimeMenu, 'String', mat2);
   set(handles.SurrMapStopTimeMenu, 'Value', 5);
   
   set(handles.GradStartTimeMenu, 'String', mat1);
   set(handles.GradStartTimeMenu, 'Value', 4);
   set(handles.GradStopTimeMenu, 'String', mat2);
   set(handles.GradStopTimeMenu, 'Value', 5);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
   handles = guihandles(gcbo);
   switch(varargin{1})
   case 'HDisp_Good'
      set(handles.DispGoodTrialsButton, 'Value', 1);
      set(handles.DispBadTrialsButton, 'Value', 0);
      
   case 'HDisp_Bad'
      set(handles.DispGoodTrialsButton, 'Value', 0);
      set(handles.DispBadTrialsButton, 'Value', 1);
      
   case 'SurrMap_Good'
      set(handles.SurrMapGoodTrialsButton, 'Value', 1);
      set(handles.SurrMapBadTrialsButton, 'Value', 0);      
      
   case 'SurrMap_Bad'
      set(handles.SurrMapGoodTrialsButton, 'Value', 0);
      set(handles.SurrMapBadTrialsButton, 'Value', 1);      
      
   case 'Grad_Good'
      set(handles.GradGoodTrialsButton, 'Value', 1);
      set(handles.GradBadTrialsButton, 'Value', 0);      
      
   case 'Grad_Bad'
      set(handles.GradGoodTrialsButton, 'Value', 0);
      set(handles.GradBadTrialsButton, 'Value', 1);      
     
   case 'HDisp_Browse'
      PATH = get(handles.PathEdit, 'String');
      Root_File = get(handles.FileNameRootEdit, 'String');
      t = sprintf('%s%s*.htb', PATH, Root_File);
      [HDisp_fname, pname] = uigetfile(t,'Select H Disparity data file');	
      if (HDisp_fname ~= 0)	
         r_ind = find(HDisp_fname == 'r');
         root = HDisp_fname(1:r_ind);
         set(handles.FileNameRootEdit, 'String', root);
         set(handles.PathEdit, 'String', pname);
         set(handles.DispEdit, 'String', HDisp_fname);
         
         %activate the surround mapping section
         set(handles.SurrMapEdit, 'Enable', 'on');
         set(handles.SurrMapBrowseButton, 'Enable', 'on');
         set(handles.SurrMapInfoButton, 'Enable', 'on');
         set(handles.SurrMapStartTimeMenu, 'Enable', 'on');
         set(handles.SurrMapStopTimeMenu, 'Enable', 'on');
         set(handles.SurrMapStartOffsetEdit, 'Enable', 'on');
         set(handles.SurrMapStopOffsetEdit, 'Enable', 'on');
         set(handles.SurrMapStartTrialsEdit, 'Enable', 'on');
         set(handles.SurrMapStopTrialsEdit, 'Enable', 'on');
         set(handles.SurrMapGoodTrialsButton, 'Enable', 'on');
         set(handles.SurrMapBadTrialsButton, 'Enable', 'on');
         set(handles.SurrMapSpikeMenu, 'Enable', 'on');
      end
      
   case 'SurrMap_Browse'
      PATH = get(handles.PathEdit, 'String');
      Root_File = get(handles.FileNameRootEdit, 'String');
      t = sprintf('%s%s*.htb', PATH, Root_File);
      [SurrMap_fname, pname] = uigetfile(t,'Select H Disparity data file');	
      if (SurrMap_fname ~= 0)	
         set(handles.PathEdit, 'String', pname);
         set(handles.SurrMapEdit, 'String', SurrMap_fname);
         
         %activate the surround mapping section
         set(handles.GradEdit, 'Enable', 'on');
         set(handles.GradBrowseButton, 'Enable', 'on');
         set(handles.GradInfoButton, 'Enable', 'on');
         set(handles.GradStartTimeMenu, 'Enable', 'on');
         set(handles.GradStopTimeMenu, 'Enable', 'on');
         set(handles.GradStartOffsetEdit, 'Enable', 'on');
         set(handles.GradStopOffsetEdit, 'Enable', 'on');
         set(handles.GradStartTrialsEdit, 'Enable', 'on');
         set(handles.GradStopTrialsEdit, 'Enable', 'on');
         set(handles.GradGoodTrialsButton, 'Enable', 'on');
         set(handles.GradBadTrialsButton, 'Enable', 'on');
         set(handles.GradSpikeMenu, 'Enable', 'on');
      end
      
   case 'Grad_Browse'
      PATH = get(handles.PathEdit, 'String');
      Root_File = get(handles.FileNameRootEdit, 'String');
      t = sprintf('%s%s*.htb', PATH, Root_File);
      [Grad_fname, pname] = uigetfile(t,'Select Disparity Gradient data file');	
      if (Grad_fname ~= 0)	
         set(handles.PathEdit, 'String', pname);
         set(handles.GradEdit, 'String', Grad_fname);
         
         %activate load data button
         set(handles.LoadButton, 'Enable', 'on');
      end      
   
   case 'H Disp file info'
      PATH = get(handles.PathEdit, 'String');
      FILE = get(handles.DispEdit, 'String');
      [return_val, info_text] = GetHTBInfo(PATH,FILE);
      if (return_val == -1)  %problem opening file
         listtext{length(listtext)+1} = 'ERROR: File could not be opened';
         set(handles.MessageListBox, 'String', listtext);
      else
         set(handles.MessageListBox, 'String', info_text);
      end
      
   case 'Surr Map file info'
      PATH = get(handles.PathEdit, 'String');
      FILE = get(handles.SurrMapEdit, 'String');
      [return_val, info_text] = GetHTBInfo(PATH,FILE);
      if (return_val == -1)  %problem opening file
         listtext{length(listtext)+1} = 'ERROR: File could not be opened';
         set(handles.MessageListBox, 'String', listtext);
      else
         set(handles.MessageListBox, 'String', info_text);
      end
      
   case 'Grad file info'
      PATH = get(handles.PathEdit, 'String');
      FILE = get(handles.GradEdit, 'String');
      [return_val, info_text] = GetHTBInfo(PATH,FILE);
      if (return_val == -1)  %problem opening file
         listtext{length(listtext)+1} = 'ERROR: File could not be opened';
         set(handles.MessageListBox, 'String', listtext);
      else
         set(handles.MessageListBox, 'String', info_text);
      end
      
   case 'load data'
      %first load hdisp data
      PATH = get(handles.PathEdit, 'String');
      FILE = get(handles.DispEdit, 'String');
      listtext = [];
      listtext{1} = 'Loading Disparity Data...';
      set(handles.MessageListBox, 'String', listtext);
      [HDisp_return_val, HDisp_good_data, HDisp_bad_data] = LoadTEMPOData(PATH,FILE);
    
      if (HDisp_return_val == -1)  %problem opening file
         listtext{length(listtext)+1} = 'ERROR: File could not be opened';
         set(handles.MessageListBox, 'String', listtext);
      else	%data are loaded, now setup some interface items
         listtext{1} = 'Loading Disparity Data...DONE';
         set(handles.MessageListBox, 'String', listtext);

         %now, set up the Start Time and Stop Time Popup Menus
         set(handles.DispStartTimeMenu, 'String', event_names);		
         set(handles.DispStartTimeMenu, 'Value', 4);		%default to Visual Stim ON
         set(handles.DispStopTimeMenu, 'String', event_names);	
         set(handles.DispStopTimeMenu, 'Value', 5);		%default to Visual Stim OFF
         %set beginning and ending trial defaults        
         set(handles.DispStartTrialsEdit, 'String', '1');		
         set(handles.DispStopTrialsEdit, 'String', num2str(size(HDisp_good_data.event_data,3)));		
         %set start and stop timing offsets for analysis - BJP 3/1/00
         set(handles.DispStartOffsetEdit, 'String', '0');
         set(handles.DispStopOffsetEdit, 'String', '0');          
         %set-up spike-channel popup menu
         good_chan_num = size(HDisp_good_data.spike_data);
         nchan = good_chan_num(1); %good_data.htb_header{SPIKE_DB}.nchannels;
         popup_text = [];
         for i=1:nchan
            spks = HDisp_good_data.spike_data(i,:,:);
            num_events = sum(sum(spks));	%count up all spikes in this spike channel
            temp = sprintf('%d   (%d events)', i, num_events);
            popup_text{length(popup_text)+1} = temp;
         end
         set(handles.DispSpikeMenu, 'String', popup_text);
      end
      
      %now load surround mapping data
      PATH = get(handles.PathEdit, 'String');
      FILE = get(handles.SurrMapEdit, 'String');
      listtext{2} = 'Loading Surround Mapping Data...';
      set(handles.MessageListBox, 'String', listtext);
      [SurrMap_return_val, SurrMap_good_data, SurrMap_bad_data] = LoadTEMPOData(PATH,FILE);
    
      if (SurrMap_return_val == -1)  %problem opening file
         listtext{length(listtext)+1} = 'ERROR: File could not be opened';
         set(handles.MessageListBox, 'String', listtext);
      else	%data are loaded, now setup some interface items
         listtext{2} = 'Loading Surround Mapping Data...DONE';
         set(handles.MessageListBox, 'String', listtext);

         %now, set up the Start Time and Stop Time Popup Menus
         set(handles.SurrMapStartTimeMenu, 'String', event_names);		
         set(handles.SurrMapStartTimeMenu, 'Value', 4);		%default to Visual Stim ON
         set(handles.SurrMapStopTimeMenu, 'String', event_names);	
         set(handles.SurrMapStopTimeMenu, 'Value', 5);		%default to Visual Stim OFF
         %set beginning and ending trial defaults        
         set(handles.SurrMapStartTrialsEdit, 'String', '1');		
         set(handles.SurrMapStopTrialsEdit, 'String', num2str(size(SurrMap_good_data.event_data,3)));		
         %set start and stop timing offsets for analysis - BJP 3/1/00
         set(handles.SurrMapStartOffsetEdit, 'String', '0');
         set(handles.SurrMapStopOffsetEdit, 'String', '0');          
         %set-up spike-channel popup menu
         good_chan_num = size(SurrMap_good_data.spike_data);
         nchan = good_chan_num(1); %good_data.htb_header{SPIKE_DB}.nchannels;
         popup_text = [];
         for i=1:nchan
            spks = SurrMap_good_data.spike_data(i,:,:);
            num_events = sum(sum(spks));	%count up all spikes in this spike channel
            temp = sprintf('%d   (%d events)', i, num_events);
            popup_text{length(popup_text)+1} = temp;
         end
         set(handles.SurrMapSpikeMenu, 'String', popup_text);
      end

      %now load gradient data
      PATH = get(handles.PathEdit, 'String');
      FILE = get(handles.GradEdit, 'String');
      listtext{3} = 'Loading Gradient Data...';
      set(handles.MessageListBox, 'String', listtext);
      [Grad_return_val, Grad_good_data, Grad_bad_data] = LoadTEMPOData(PATH,FILE);
    
      if (Grad_return_val == -1)  %problem opening file
         listtext{length(listtext)+1} = 'ERROR: File could not be opened';
         set(handles.MessageListBox, 'String', listtext);
      else	%data are loaded, now setup some interface items
         listtext{3} = 'Loading Gradient Data...DONE';
         set(handles.MessageListBox, 'String', listtext);

         %now, set up the Start Time and Stop Time Popup Menus
         set(handles.GradStartTimeMenu, 'String', event_names);		
         set(handles.GradStartTimeMenu, 'Value', 4);		%default to Visual Stim ON
         set(handles.GradStopTimeMenu, 'String', event_names);	
         set(handles.GradStopTimeMenu, 'Value', 5);		%default to Visual Stim OFF
         %set beginning and ending trial defaults        
         set(handles.GradStartTrialsEdit, 'String', '1');		
         set(handles.GradStopTrialsEdit, 'String', num2str(size(Grad_good_data.event_data,3)));		
         %set start and stop timing offsets for analysis - BJP 3/1/00
         set(handles.GradStartOffsetEdit, 'String', '0');
         set(handles.GradStopOffsetEdit, 'String', '0');          
         %set-up spike-channel popup menu
         good_chan_num = size(Grad_good_data.spike_data);
         nchan = good_chan_num(1); %good_data.htb_header{SPIKE_DB}.nchannels;
         popup_text = [];
         for i=1:nchan
            spks = Grad_good_data.spike_data(i,:,:);
            num_events = sum(sum(spks));	%count up all spikes in this spike channel
            temp = sprintf('%d   (%d events)', i, num_events);
            popup_text{length(popup_text)+1} = temp;
         end
         set(handles.GradSpikeMenu, 'String', popup_text);
      end      
      
      if((Grad_return_val == -1) | (SurrMap_return_val == -1) | (HDisp_return_val == -1))
         listtext{length(listtext)+1} = 'ERROR: One File could not be opened, abort analysis';
         set(handles.MessageListBox, 'String', listtext);
      else
         set(handles.AnalyzeButton, 'Enable', 'on');
      end
      
   case 'analyze data'
      DispStartCode = get(handles.DispStartTimeMenu, 'Value');
      DispStopCode = get(handles.DispStopTimeMenu, 'Value');
      DispBegTrial = str2num(get(handles.DispStartTrialsEdit, 'String'));
      DispEndTrial = str2num(get(handles.DispStopTrialsEdit, 'String'));
      DispStartOffset = str2num(get(handles.DispStartOffsetEdit, 'String'));
      DispStopOffset = str2num(get(handles.DispStopOffsetEdit, 'String'));
    
      DispSpikeChan = get(handles.DispSpikeMenu, 'Value');
      PATH = get(handles.PathEdit, 'String');
      DispFILE = get(handles.DispEdit, 'String');
      %depending on how the Good Trials/Bad Trials radio buttons are set, analyze the appropriate set of trials
      good_flag = get(handles.DispGoodTrialsButton, 'Value');
      if (good_flag)
          Disp_select_data = HDisp_good_data;
      else
          Disp_select_data = HDisp_bad_data;
      end

      SurrMapStartCode = get(handles.SurrMapStartTimeMenu, 'Value');
      SurrMapStopCode = get(handles.SurrMapStopTimeMenu, 'Value');
      SurrMapBegTrial = str2num(get(handles.SurrMapStartTrialsEdit, 'String'));
      SurrMapEndTrial = str2num(get(handles.SurrMapStopTrialsEdit, 'String'));
      SurrMapStartOffset = str2num(get(handles.SurrMapStartOffsetEdit, 'String'));
      SurrMapStopOffset = str2num(get(handles.SurrMapStopOffsetEdit, 'String'));
    
      SurrMapSpikeChan = get(handles.SurrMapSpikeMenu, 'Value');
      PATH = get(handles.PathEdit, 'String');
      SurrMapFILE = get(handles.SurrMapEdit, 'String');
      %depending on how the Good Trials/Bad Trials radio buttons are set, analyze the appropriate set of trials
      good_flag = get(handles.SurrMapGoodTrialsButton, 'Value');
      if (good_flag)
          SurrMap_select_data = SurrMap_good_data;
      else
          SurrMap_select_data = SurrMap_bad_data;
      end

      GradStartCode = get(handles.GradStartTimeMenu, 'Value');
      GradStopCode = get(handles.GradStopTimeMenu, 'Value');
      GradBegTrial = str2num(get(handles.GradStartTrialsEdit, 'String'));
      GradEndTrial = str2num(get(handles.GradStopTrialsEdit, 'String'));
      GradStartOffset = str2num(get(handles.GradStartOffsetEdit, 'String'));
      GradStopOffset = str2num(get(handles.GradStopOffsetEdit, 'String'));
    
      GradSpikeChan = get(handles.GradSpikeMenu, 'Value');
      PATH = get(handles.PathEdit, 'String');
      GradFILE = get(handles.GradEdit, 'String');
      %depending on how the Good Trials/Bad Trials radio buttons are set, analyze the appropriate set of trials
      good_flag = get(handles.GradGoodTrialsButton, 'Value');
      if (good_flag)
          Grad_select_data = Grad_good_data;
      else
          Grad_select_data = Grad_bad_data;
      end      
      
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
end %end switch


end
