function varargout = roi_crawl(varargin)
% GUI to view ROI plots
%       usage:
%       roi_crawl(basename, intensity_file)
%
%       example:
%       roi_crawl('2006-02-02-2_reg.imagine', '2006-02-02-2_reg_imagine_roi_intensities.intensity');
%
%   this script uses roi_traceviewer script. the options to which have to be
%   changed in this code directly. Have to add functionality to do that from
%   the command line. 
%   function roi_tracemaker at the end of the script is same as
%   roi_traceviewer. Its defaults can be changed to change things such as
%   filtering, plotting parameters. The help is at the bottom neat the
%   function
%
% Diwakar Turaga 2006-10-15 (my first GUI)
%
% ROI_CRAWL M-file for roi_crawl.fig
%      ROI_CRAWL, by itself, creates a new ROI_CRAWL or raises the existing
%      singleton*.
%
%      H = ROI_CRAWL returns the handle to a new ROI_CRAWL or the handle to
%      the existing singleton*.
%
%      ROI_CRAWL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROI_CRAWL.M with the given input arguments.
%
%      ROI_CRAWL('Property','Value',...) creates a new ROI_CRAWL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before roi_crawl_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to roi_crawl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help roi_crawl

% Last Modified by GUIDE v2.5 15-Oct-2006 07:17:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @roi_crawl_OpeningFcn, ...
                   'gui_OutputFcn',  @roi_crawl_OutputFcn, ...
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


% --- Executes just before roi_crawl is made visible.
function roi_crawl_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to roi_crawl (see VARARGIN)

% set_new_roi_number(roi_number);

% Choose default command line output for roi_crawl
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes roi_crawl wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% basename = '2006-02-02-2_reg_best_2sided.imagine';
% intensity_file = '2006-02-02-2_reg_best_2sided.imagine_roi_intensities.intensity';

if ~isempty(varargin)
    basename = varargin{1};
    intensity_file = varargin{2};
% else  
%     basename = '2006-02-02-2_reg_best_2sided.imagine';
%     intensity_file = '2006-02-02-2_reg_best_2sided.imagine_roi_intensities.intensity';
end

options.filterpaneldata = 0;
[intensities_filt, inten_normal_unfilt, stimulus_unfilt, labels, stimrange, timewindowpretrans] = roi_tracemaker(basename, intensity_file, options);
options.filterpaneldata = 1;
[intensities_filt, inten_normal_unfilt, stimulus_filt, labels, stimrange, timewindowpretrans] = roi_tracemaker(basename, intensity_file, options);

handles.intensities_filt = intensities_filt;
handles.inten_normal_unfilt = inten_normal_unfilt;
handles.stimulus_filt = stimulus_filt;
handles.stimulus_unfilt = stimulus_unfilt;
handles.labels = labels;
handles.stimrange = stimrange;
handles.timewindowpretrans = timewindowpretrans;

handles.stimulus = handles.stimulus_unfilt; % default
roi_number = 1;
roi_number = plot_roi(roi_number, handles);
set(handles.edit1,'String', num2str(roi_number))

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = roi_crawl_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roi_number = str2double(get(handles.edit1, 'String'));
roi_number = roi_number + 1;
roi_number = plot_roi(roi_number, handles);
set(handles.edit1,'String', num2str(roi_number))

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roi_number = str2double(get(handles.edit1, 'String'));
roi_number = roi_number - 1;
roi_number = plot_roi(roi_number, handles);
set(handles.edit1,'String', num2str(roi_number))

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

roi_number = str2double(get(hObject, 'String'));
roi_number = plot_roi(roi_number, handles);
set(handles.edit1,'String', num2str(roi_number))

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    handles.stimulus = handles.stimulus_filt;
    roi_number = str2double(get(handles.edit1, 'String'));
    roi_number = plot_roi(roi_number, handles);
    
elseif button_state == get(hObject,'Min')
    handles.stimulus = handles.stimulus_unfilt;
    roi_number = str2double(get(handles.edit1, 'String'));
    roi_number = plot_roi(roi_number, handles);
end

%--------------------------------------------
function roi_number = plot_roi(roi_number, handles)

button_state = get(handles.checkbox1,'Value');
if button_state == get(handles.checkbox1,'Max')
    handles.stimulus = handles.stimulus_filt; 
elseif button_state == get(handles.checkbox1,'Min')
    handles.stimulus = handles.stimulus_unfilt;
end

intensities_filt = handles.intensities_filt;
inten_normal_unfilt = handles.inten_normal_unfilt;
stimulus = handles.stimulus;
labels = handles.labels;
stimrange = handles.stimrange;
timewindowpretrans = handles.timewindowpretrans;

roi_number = round(roi_number);

if roi_number > size(intensities_filt, 1)
    roi_number = size(intensities_filt, 1);
end

if roi_number < 1
    roi_number = 1;
end

axes(handles.filtered)
plot(intensities_filt(roi_number, :),'color',[0 0 0])

axes(handles.unfiltered)
plot(inten_normal_unfilt(roi_number, :))

trialmap = trialcolor(size(stimulus,1));
if 1 <= size(stimulus, 4)
    axes(handles.stim1)
    for i = 1:size(stimulus, 1)
        plot(squeeze(stimulus(i, roi_number, :, 1)),'color',trialmap{i});
        axis tight; hold on; title(labels{1});
        hold on
    end
    meanstim = (squeeze(nanmean(stimulus(:,roi_number,:,1),1)))';
    plot(meanstim,'color',[0 0 0],'LineWidth',2)
    yaxrange = get(gca,'Ylim');
    xaxrange = get(gca,'Xlim');
    tracehand = get(gca,'Children');
    hrect = rectangle('Position',[timewindowpretrans+1, yaxrange(1), stimrange, yaxrange(2)-yaxrange(1)],'EdgeColor','none','FaceColor',[0.97 0.97 0.97]);
    set(gca,'Xlim',xaxrange,'Ylim',yaxrange);
    set(gca,'Children',[tracehand;hrect']);
    hold off
end

% first stimulus is usually KCl,

y_lim = [min(min((squeeze(min(stimulus(:, roi_number, :, 2:end))))')) max(max((squeeze(max(stimulus(:, roi_number, :, 2:end))))'))];

if 2 <= size(stimulus, 4)
    axes(handles.stim2)
    for i = 1:size(stimulus, 1)
        plot(squeeze(stimulus(i, roi_number, :, 2)),'color',trialmap{i}); 
        axis tight; hold on; title(labels{2})
    end
    meanstim = (squeeze(nanmean(stimulus(:,roi_number,:,2),1)))';
    plot(meanstim,'color',[0 0 0],'LineWidth',2)
    ylim(y_lim); 
    yaxrange = get(gca,'Ylim');
    xaxrange = get(gca,'Xlim');
    tracehand = get(gca,'Children');
    hrect = rectangle('Position',[timewindowpretrans+1, yaxrange(1), stimrange, yaxrange(2)-yaxrange(1)],'EdgeColor','none','FaceColor',[0.97 0.97 0.97]);
    set(gca,'Xlim',xaxrange,'Ylim',yaxrange);
    set(gca,'Children',[tracehand;hrect']);
    hold off
end

if 3 <= size(stimulus, 4)
    axes(handles.stim3)
    for i = 1:size(stimulus, 1)
        plot(squeeze(stimulus(i, roi_number, :, 3)),'color',trialmap{i}); 
        ylim(y_lim); axis tight; hold on; title(labels{3});
    end
    meanstim = (squeeze(nanmean(stimulus(:,roi_number,:,3),1)))';
    plot(meanstim,'color',[0 0 0],'LineWidth',2)
    ylim(y_lim); 
    yaxrange = get(gca,'Ylim');
    xaxrange = get(gca,'Xlim');
    tracehand = get(gca,'Children');
    hrect = rectangle('Position',[timewindowpretrans+1, yaxrange(1), stimrange, yaxrange(2)-yaxrange(1)],'EdgeColor','none','FaceColor',[0.97 0.97 0.97]);
    set(gca,'Xlim',xaxrange,'Ylim',yaxrange);
    set(gca,'Children',[tracehand;hrect']);
    hold off
end

if 4 <= size(stimulus, 4)
    axes(handles.stim4)
    for i = 1:size(stimulus, 1)
        plot(squeeze(stimulus(i, roi_number, :, 4)),'color',trialmap{i});
        ylim(y_lim); axis tight; hold on; title(labels{4});
    end
    meanstim = (squeeze(nanmean(stimulus(:,roi_number,:,4),1)))';
    plot(meanstim,'color',[0 0 0],'LineWidth',2)
    ylim(y_lim); 
    yaxrange = get(gca,'Ylim');
    xaxrange = get(gca,'Xlim');
    tracehand = get(gca,'Children');
    hrect = rectangle('Position',[timewindowpretrans+1, yaxrange(1), stimrange, yaxrange(2)-yaxrange(1)],'EdgeColor','none','FaceColor',[0.97 0.97 0.97]);
    set(gca,'Xlim',xaxrange,'Ylim',yaxrange);
    set(gca,'Children',[tracehand;hrect']);
    hold off
end

if 5 <= size(stimulus, 4)
    axes(handles.stim5)
    for i = 1:size(stimulus, 1)
        plot(squeeze(stimulus(i, roi_number, :, 5)),'color',trialmap{i});
        axis tight; hold on; title(labels{5});
    end
    meanstim = (squeeze(nanmean(stimulus(:,roi_number,:,5),1)))';
    plot(meanstim,'color',[0 0 0],'LineWidth',2)
    ylim(y_lim); 
    yaxrange = get(gca,'Ylim');
    xaxrange = get(gca,'Xlim');
    tracehand = get(gca,'Children');
    hrect = rectangle('Position',[timewindowpretrans+1, yaxrange(1), stimrange, yaxrange(2)-yaxrange(1)],'EdgeColor','none','FaceColor',[0.97 0.97 0.97]);
    set(gca,'Xlim',xaxrange,'Ylim',yaxrange);
    set(gca,'Children',[tracehand;hrect']);
    hold off
end

if 6 <= size(stimulus, 4)
    axes(handles.stim6)
    for i = 1:size(stimulus, 1)
        plot(squeeze(stimulus(i, roi_number, :, 6)),'color',trialmap{i});
        axis tight; hold on; title(labels{6})
        hold on
    end
    meanstim = (squeeze(nanmean(stimulus(:,roi_number,:,6),1)))';
    plot(meanstim,'color',[0 0 0],'LineWidth',2)
    ylim(y_lim); 
    yaxrange = get(gca,'Ylim');
    xaxrange = get(gca,'Xlim');
    tracehand = get(gca,'Children');
    hrect = rectangle('Position',[timewindowpretrans+1, yaxrange(1), stimrange, yaxrange(2)-yaxrange(1)],'EdgeColor','none','FaceColor',[0.97 0.97 0.97]);
    set(gca,'Xlim',xaxrange,'Ylim',yaxrange);
    set(gca,'Children',[tracehand;hrect']);
    hold off
end

if 7 <= size(stimulus, 4)
    axes(handles.stim7)
    for i = 1:size(stimulus, 1)
        plot(squeeze(stimulus(i, roi_number, :, 7)),'color',trialmap{i});
        axis tight; hold on; title(labels{7})
        hold on
    end
    meanstim = (squeeze(nanmean(stimulus(:,roi_number,:,7),1)))';
    plot(meanstim,'color',[0 0 0],'LineWidth',2)
    ylim(y_lim); 
    yaxrange = get(gca,'Ylim');
    xaxrange = get(gca,'Xlim');
    tracehand = get(gca,'Children');
    hrect = rectangle('Position',[timewindowpretrans+1, yaxrange(1), stimrange, yaxrange(2)-yaxrange(1)],'EdgeColor','none','FaceColor',[0.97 0.97 0.97]);
    set(gca,'Xlim',xaxrange,'Ylim',yaxrange);
    set(gca,'Children',[tracehand;hrect']);
    hold off
end

if 8 <= size(stimulus, 4)
    axes(handles.stim8)
    for i = 1:size(stimulus, 1)
        plot(squeeze(stimulus(i, roi_number, :, 8)),'color',trialmap{i});
        axis tight; hold on; title(labels{8})
        hold on
    end
    meanstim = (squeeze(nanmean(stimulus(:,roi_number,:,8),1)))';
    plot(meanstim,'color',[0 0 0],'LineWidth',2)
    ylim(y_lim); 
    yaxrange = get(gca,'Ylim');
    xaxrange = get(gca,'Xlim');
    tracehand = get(gca,'Children');
    hrect = rectangle('Position',[timewindowpretrans+1, yaxrange(1), stimrange, yaxrange(2)-yaxrange(1)],'EdgeColor','none','FaceColor',[0.97 0.97 0.97]);
    set(gca,'Xlim',xaxrange,'Ylim',yaxrange);
    set(gca,'Children',[tracehand;hrect']);
    hold off
end




%----------------------------------------------


function [intensities_filt, inten_normal_unfilt, stimulus, labels, stimrange, timewindowpretrans] = roi_tracemaker(basename, intensityfilename, options)
%roi_traceviewer takes the roi intensity matrix calculated by
%calc_stack_roi_inten and returns a set of panels with useful plots. As the
%user hits any key the panels refresh to show information from the next
%sequential ROI. The raw data, a bessel filtered version of the data, and
%trial overlays for each stimulus (plotted using trialcolor for
%color picking)as well as the mean response for all trials.
% Syntax:
%   roi_traceviewer(basename, intensityfilename)
%   roi_traceviewer(basename, intensityfilename, options)
% where
%   basename is the original basic filename from the imagine experimental data (also the header name), 
%             BTW:  the header file from the original data should automatically be saved
%             as a .mat or .imagine file into the parent directory for this job by
%             the previous processing function (registration, etc.). If not, you should do this.
%             NOTE: if the previous processing function modified the shape
%             of the data (e.g. dropped the first and last stacks) the
%             copied header needs to be edited to reflect this (you may need to check that
%             this is so).  
%   intensityfilename is the name of the saved intensity matrix file created by
%             calc_stack_roi_inten 
%   options is a structure indicating some of the extra behavior you would
%             or would not like roi_traceviewer to exhibit.
%        .timewindowpretrans can be set to change the pre-stimulus period
%             used for the trial snippeting time window.
%             Default is 4 (stacks).
%        .timewindowposttrans can be set to change the after stimulus-transition period
%             used for the trial snippeting time window.
%             Default is 15 (stacks).
%        .backgrndmeanlength can be set to change the pre-stimulus background period
%             used for normalization.
%             Default is 3 (stacks).
%        .n is the number of poles to used for the bessel noise reduction
%             filter. Default is 3.
%        .Wn Wn is a two-element vector, Wn = [w1 w2], butter returns an order 2*n 
%             digital bandpass filter with passband w1 < w < w2. Default is [0.02 0.5].
%        .filt_type indicates the type of filter you want to use. Options
%             are: 'high', 'low', 'stop'.  If you say nothing, the automatic
%             default is bandpass.
%        .filterpaneldata indicates whether you want the filtered or
%             unfiltered data presented in the trial overlay panels.
%             Default is 1(true) for presenting the filtered data.
%Copywrite 2006 by Diwakar Turaga & Terrence Holekamp

if nargin < 3
    options.timewindowpretrans = 4;
    options.timewindowposttrans = 15;
    options.backgrndmeanlength = 3;
    options.n = 5;
    options.Wn = [0.45];
    options.paneldatatype = 1;
end

if(~isfield(options, 'timewindowpretrans'))
    options.timewindowpretrans = 4;
end
if(~isfield(options, 'timewindowposttrans'))
    options.timewindowposttrans = 15;
end
if(~isfield(options, 'backgrndmeanlength'))
    options.backgrndmeanlength = 3;
end
if(~isfield(options, 'n'))
    options.n = 5;
end
if(~isfield(options, 'Wn'))
    options.Wn = [0.45];
end
if(~isfield(options, 'filterpaneldata'))
    options.filterpaneldata = 1;
end

timewindowpretrans = options.timewindowpretrans;
timewindowposttrans = options.timewindowposttrans;
bkglng = options.backgrndmeanlength;
n = options.n;
Wn = options.Wn;
filtpanel = options.filterpaneldata;

%[b,a] = besself(n,Wo);
if(isfield(options, 'filt_type'))
    [b,a] = butter(n,Wn,options.filt_type);
else
    [b,a] = butter(n,Wn,'low');
end

if ~exist('intensities','var')
    load (intensityfilename,'-mat')
end

smm = stackmm(basename);
header = smm.header;
stimuli = header.stim_lookup;
labels = header.stim_labels;

size_inten = size(intensities);
number_of_ROI  = size_inten(1);
number_of_stacks = size_inten(2);
% number_of_trials = max(unique(header.trial_lookup(~isnan(header.trial_lookup))));
number_of_diff_stimuli = max(unique(stimuli));

number_of_trials = 0; % this is so that if we dont register the whole stacks, we can find out exactly how many of the trials had been registered
for i = 1:number_of_stacks
    if stimuli(i) == 1 && stimuli(i-1) == 0
        number_of_trials = number_of_trials + 1;
    end
end

% filter the data
intensities(:,1) = intensities(:,2);
intensities_filt = filtfilt(b,a,intensities');
intensities_filt = intensities_filt';

for i = 1:number_of_ROI
    inten_avg_unfilt = nanmean(intensities(i,:));
    inten_normal_unfilt(i,:) = (intensities(i,:) - inten_avg_unfilt)/inten_avg_unfilt;
end

if bkglng > timewindowpretrans
    disp('you may be asking for too long a background averaging period')
end

if filtpanel
    data = intensities_filt;
else
    data = intensities;
end

stimulus = NaN(number_of_trials,number_of_ROI, timewindowpretrans+timewindowposttrans+1, number_of_diff_stimuli);

for i = 1:number_of_diff_stimuli
    trial_counter(i) = 1;
end

i = 1;
 
while i <=number_of_stacks-1
    
    if (stimuli(i)==0 && stimuli(i+1)~=0)
            
        for j = 1:number_of_diff_stimuli
           
            if stimuli(i+1) == j
                for k = 1:number_of_ROI
                    inten_backgrnd_avg = nanmean(data(k,i-bkglng:i));
                    if i+timewindowposttrans <= size(data,2)
                        inten_snip = data(k,i-timewindowpretrans:i+timewindowposttrans);
                    else
                         inten_snip = data(k,i-timewindowpretrans:end);
                    end
                    inten_normal = (inten_snip - inten_backgrnd_avg)/inten_backgrnd_avg;
                    stimulus(trial_counter(j),k,1:length(inten_normal),j) = inten_normal;
                end
                i = i+1;
                trial_counter(j) = trial_counter(j)+1;
            end 
            
        end
                        
            
    else
        i=i+1;
    end
end

r = diff(find(header.stim_lookup ==1));
r(r == 1) = 0;
s = find(r);
stimrange = s(1);



