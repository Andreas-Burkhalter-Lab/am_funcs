function merec_gui_thresholding(options)
% MEREC_GUI_THRESHOLDING: a gui-based, wysiwyg tool assisting choices 
%                         of bandpass filtering, conditioning filters,
%                         detection filters, and thresholds before 
%                         snippeting 
%
%      use: merec_gui_thresholding;
%           Run from anywhere, preferably the directory where files
%           to be analyzed are stored
%      
%      option values: not implemented
%
%      key features: merec_gui_thresholding is an alternate way to scan
%                    data and manually choose preferred snippeting options
%
%
%      SEE ALSO: SNIPPETFILE, SNIPOPTIONS, MEREC2VLV 
%
%      Copyright Julian P. Meeks, 2007
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
% Date created: 2007_02_07
% Last modified: 2007_05_14
% Last editor: JPM
% To-do list: ---------
% -----------------------------------------------------------------------

% -- Initialization goes here ---

% -- Initialization end ---------

version_id = 1.0;

% Initialize masterfig, hiding for now --------------------------------
masterfig = figure('Visible', 'on', 'Name', 'MEREC Thresholding', ...
                   'Position', [0 24 1200 700]); 
% ---------------------------------------------------------------------

% Initialize orig_axis to contain original data input -----------------
orig_panel = uipanel(masterfig, 'Title', 'Original Waveform', ...
               'Units', 'pixels', 'Position', [10,310,420,280], ...
               'BackgroundColor', 'white');
orig_axis = axes('Parent', orig_panel, 'Units', 'pixels', ...
               'Position', [36,36,360,212]);
orig_xlabel = uicontrol(orig_panel, 'Style', 'text',...
                        'String', 'seconds', 'Units', 'normalized',...
                        'Position', [0.5,0.01,.1,.05], 'BackgroundColor',...
                        'white');
orig_zoombtn = uicontrol(orig_panel, 'String', 'Zoom', ...
               'Position', [420,246,50,20], ...
               'Callback', {@zoombtn_Callback});
orig_panbtn = uicontrol(orig_panel, 'String', 'Pan', ...
               'Position', [420,222,50,20], ...
               'Callback', {@panbtn_Callback});
orig_envbtn = uicontrol(orig_panel, 'String', 'Envelope', ...
               'Position', [420,198,60,20], ...
               'Callback', {@envbtn_Callback});
orig_resetbtn = uicontrol(orig_panel, 'String', 'Reset', ...
               'Position', [420,36,60,20], ...
               'Callback', {@resetbtn_Callback});
% ---------------------------------------------------------------------

% Initialize mod_axis to contain modified data ------------------------
mod_panel = uipanel(masterfig, 'Title', 'Modified Waveform', ...
               'Units', 'pixels', 'Position', [10,10,420,290], ...
               'BackgroundColor', 'white');
mod_axis = axes('Parent', mod_panel, 'Units', 'pixels', ...
               'Position', [36,36,360,212]);
mod_xlabel = uicontrol(mod_panel, 'Style', 'text',...
                        'String', 'seconds', 'Units', 'normalized',...
                        'Position', [0.5,0.01,.1,.05], 'BackgroundColor',...
                        'white');
mod_threshline = 0;

% ---------------------------------------------------------------------

% Initialize bandpass_panel panel to contain bandpass info ------------
bandpass_panel = uipanel(masterfig, 'Title', 'Bandpass Filtering', ...
               'Units', 'pixels','BackgroundColor', 'cyan', ...
               'Position', [500,470,290,120]);
bandpass_check = uicontrol(bandpass_panel, 'Style', 'checkbox',...
               'Units', 'Pixels', 'Position', [-8,104,16,16], 'Value', 0,...
               'Min', 0, 'Max', 1);
bandpass_prompt = uicontrol(bandpass_panel, 'Style', 'text',...
               'Units', 'Pixels', 'Position', [2,80,280,16],...
               'Value', 0, 'BackgroundColor', 'cyan', 'String',...
               'Please enter the limits of the bandpass filter (in Hz):');
bandpass_low = uicontrol(bandpass_panel, 'Style', 'edit', ...
               'Units', 'Pixels', 'Position', [24,60,60,16],...
               'BackgroundColor', 'white',...
               'String', 'Highpass'); % use str2num( ) to convert back to number
bandpass_wlabel = uicontrol(bandpass_panel, 'Style', 'text',...
               'BackgroundColor', 'cyan', 'String', '  <  w <  ',...
               'Units', 'Pixels', 'Position', [90,60,60,16]);
bandpass_high = uicontrol(bandpass_panel, 'Style', 'edit', ...
               'Units', 'Pixels', 'Position', [160,60,60,16],...
               'BackgroundColor', 'white',...
               'String', 'Lowpass'); % use str2num( ) to convert back to number
bandpass_poles = uibuttongroup('Parent', bandpass_panel, ...
               'Title', 'Filter Pole #','Units', 'Pixels',...
               'Position', [12,10,260,40], 'BackgroundColor','cyan');
bandpass_2 = uicontrol(bandpass_poles, 'Style', 'radiobutton',...
               'String', ' 2 ', 'Units','pixels','Position', [12,10,40,16],...
               'BackgroundColor', 'cyan');
bandpass_4 = uicontrol(bandpass_poles, 'Style', 'radiobutton',...
               'String', ' 4 ', 'Units', 'pixels', 'Position', [72,10,40,16],...
               'BackgroundColor', 'cyan');
bandpass_8 = uicontrol(bandpass_poles, 'Style', 'radiobutton',...
               'String', ' 8 ', 'Units', 'pixels', 'Position', [120,10,40,16],...
               'BackgroundColor', 'cyan');
bandpass_16 = uicontrol(bandpass_poles, 'Style', 'radiobutton',...
               'String', ' 16 ', 'Units', 'pixels', 'Position', [172,10,40,16],...
               'BackgroundColor', 'cyan');
align([bandpass_2 bandpass_4 bandpass_8 bandpass_16], 'Distribute', 'Center');
bandpass_apply = uicontrol(bandpass_panel, 'Units', 'pixels',...
                          'Position', [230,60,50,20], 'String', 'Apply',...
                          'Callback', {@bandpass_apply_Callback});
    
% ---------------------------------------------------------------------

% Initialize detection_panel panel to contain detection filter info ---
detection_panel = uipanel(masterfig, 'Title', 'Detection filter', ...
               'Units', 'pixels', 'BackgroundColor', 'green', ...
               'Position', [500, 400, 290, 60]);
det_panel_label = uicontrol(detection_panel, 'Units', 'normalized',...
               'Style', 'text', 'String', 'not yet implemented', ...
               'Position', [0.3,0.5,0.4,0.3], 'BackgroundColor', 'green');
% ---------------------------------------------------------------------

% Initialize thresholding_panel to contain thresholding info ----------
thresholding_panel = uipanel(masterfig, 'Title', 'Thresholding', ...
               'Units', 'pixels', 'BackgroundColor', 'yellow', ...
               'Position', [500,170,290,220]);
thresholding_slider = uicontrol(thresholding_panel, 'Style', 'slider',...
               'Units', 'normalized', 'Position', [0.1,0.8,0.8,0.1],...
               'Min', 2, 'Max', 12, 'Value', 6, 'SliderStep', [1/100 1/20],...
               'Callback', {@thresh_slider_move});
thresh_slider_label_def = sprintf('The current threshold scaling value is: %0.2f',...
                                 get(thresholding_slider, 'Value'));
thresholding_slider_label = uicontrol(thresholding_panel, 'Style', 'text',...
               'Units', 'normalized', 'Position', [0.1,0.9,0.8,0.08],...
               'BackgroundColor', 'yellow', 'String', thresh_slider_label_def);
thresh_mean_noise = 0;
thresh_mean_noise_def = sprintf('Mean noise: %0.5f', thresh_mean_noise);
thresh_mean_noise_label = uicontrol(thresholding_panel, 'Style', 'text', ...
               'Units', 'normalized', 'Position', [0.01,0.6,0.25,0.16],...
               'BackgroundColor', 'yellow', 'String', thresh_mean_noise_def);
thresh_value_def = sprintf('Threshold value: %0.5f', ...
                            thresh_mean_noise*get(thresholding_slider,'Value'));
thresh_value_label = uicontrol(thresholding_panel, 'Style', 'text',...
               'Units', 'normalized', 'Position', [0.26,0.6,0.3,0.16],...
               'BackgroundColor', 'yellow', 'String', thresh_value_def);
thresh_update = uicontrol(thresholding_panel, 'Units', 'normalized',...
               'Position', [0.79,0.62,0.2,0.12],'String','Update',...
               'Callback', {@find_mna});
thresh_psnip_def = sprintf('Pseudo-\nsnips: %d', 0.0);
threshold_pseudosnip_num = uicontrol(thresholding_panel, 'Units','Normalized',...
                        'Position',[0.56,0.6,0.22,0.16], 'Style','text',...
                        'BackgroundColor', 'red', 'String', thresh_psnip_def);
psnip_plot = 0;
snip_plot = 0;
thresh_divider = uicontrol(thresholding_panel, 'Units', 'normalized', ...
                        'Position', [0.01,0.58,0.98,0.04], 'Style', 'text',...
                        'BackgroundColor', 'yellow', ...
                        'String', '=======================================');
thresh_rebnd_check = uicontrol(thresholding_panel, 'Units', 'normalized',...
                        'Position', [0.01,0.48,0.08,0.08], 'Style', 'checkbox',...
                        'Value', 0, 'BackgroundColor', 'yellow');
thresh_rebnd_label = uicontrol(thresholding_panel, 'Units', 'normalized',...
                        'Position', [0.09,0.45,0.16,0.1], 'Style', 'text',...
                        'String', 'Rebound', 'BackgroundColor', 'yellow');
thresh_rebnd_S_box = uicontrol(thresholding_panel, 'Style', 'edit',...
                        'Units', 'normalized', 'Position', [0.28,0.46,0.3,0.1],...
                        'String', 'Rebound scan');
thresh_rebnd_C_box = uicontrol(thresholding_panel, 'Style', 'edit',...
                        'Units', 'normalized', 'Position', [0.6, 0.46,0.36,0.1],...
                        'String', 'Rebound curvature');
thresh_peaktrough_label = uicontrol(thresholding_panel, 'Style', 'text', ...
                        'Units', 'normalized', 'Position', [0.01,0.26,0.4,0.16],...
                        'String', '# samples between peak & trough:',...
                        'BackgroundColor', 'yellow');
thresh_peaktrough_box = uicontrol(thresholding_panel, 'Style', 'edit', ...
                        'Units', 'normalized', 'Position', [0.4,0.30,0.1,0.12],...
                        'String', '30');
thresh_close_label = uicontrol(thresholding_panel, 'Style', 'text',...
                        'Units', 'normalized', 'Position', [0.51,0.26,0.4,0.16],...
                        'String', 'min # samples between peaks:',...
                        'BackgroundColor', 'yellow');
thresh_close_box = uicontrol(thresholding_panel, 'Style', 'edit', ...
                        'Units', 'normalized', 'Position', [0.88,0.30,0.1,0.12],...
                        'String', '30');
thresh_troughdepth_label = uicontrol(thresholding_panel, 'Style', 'text',...
                        'Units', 'normalized', 'Position', [0.02,0.12,0.3,0.12],...
                        'String', 'Trough depth:', 'BackgroundColor', 'yellow');
thresh_troughdepth_box = uicontrol(thresholding_panel, 'Style', 'edit',...
                        'Units', 'normalized', 'Position', [0.3,0.16,0.12,0.12],...
                        'String', '0.5');
thresh_calcsnip_button = uicontrol(thresholding_panel, 'Units', 'normalized',...
                        'Position', [0.04,0.02,0.3,0.12], 'String', 'Calc Snips',...
                        'Callback', {@calc_snips});
thresh_snip_value = sprintf('The number of \nreal snips is: %d', 0);
thresh_snip_label = uicontrol(thresholding_panel, 'Units', 'normalized', ...
                        'Style', 'text', 'Position', [0.65,0.16,0.32,0.12],...
                        'BackgroundColor', 'red','String', thresh_snip_value);
                    
polarity_radio = uibuttongroup('Parent', thresholding_panel, ...
               'Title', 'Spike Polarity:','Units', 'normalized',...
               'Position', [0.4 0.01 0.58 0.12], 'BackgroundColor','red');
polarity_pos = uicontrol(polarity_radio, 'Style', 'radiobutton',...
               'String', '+', 'Units','normalized','Position', [0.05 0.01 0.25 0.99],...
               'BackgroundColor', 'red');
polarity_both = uicontrol(polarity_radio, 'Style', 'radiobutton',...
               'String', '+/-', 'Units','normalized','Position', [0.38 0.01 0.25 0.99],...
               'BackgroundColor', 'red', 'value', 1);
polarity_neg = uicontrol(polarity_radio, 'Style', 'radiobutton',...
               'String', '-', 'Units','normalized','Position', [0.71 0.01 0.25 0.99],...
               'BackgroundColor', 'red');

% ---------------------------------------------------------------------

% Initialize io_panel to contain basic input/output options -----------
io_panel = uipanel(masterfig, 'Title', 'Input/Output',...
               'Units', 'pixels', 'BackgroundColor', 'magenta', ...
               'Position', [500,10,290,90]);
io_write_options_button = uicontrol(io_panel, 'Units', 'normalized',...
                        'Position', [0.05,0.65,0.4,0.24], 'String', 'Write Options',...
                        'Callback', {@write_ops});
io_snip_button = uicontrol(io_panel, 'Units', 'normalized',...
                        'Position', [0.05,0.15,0.4,0.24], 'String', '.ssnp It',...
                        'Callback', {@snip_it});
io_valve_button = uicontrol(io_panel, 'Units', 'normalized', ...
                        'Position', [0.5,0.65,0.4,0.24], 'String', '.vlv It',...
                        'Callback', {@gen_vlv_file});
io_doall_button = uicontrol(io_panel, 'Units', 'normalized', ...
                        'Position', [0.5,0.15,0.4,0.24], 'String', 'Do All',...
                        'Callback', {@snip_all_in_dir});
% ---------------------------------------------------------------------

% Set UI values to normalized to maintain spatial relationships -------
set([masterfig,orig_panel,orig_axis,orig_zoombtn,orig_panbtn,...
               orig_resetbtn,orig_envbtn,...
               mod_panel,mod_axis,...
               bandpass_panel,bandpass_check, bandpass_prompt,...
               bandpass_wlabel, bandpass_low, bandpass_high,... 
               bandpass_poles, bandpass_2, bandpass_4, bandpass_8,...
               bandpass_16, bandpass_apply,...
               detection_panel,...
               thresholding_panel,io_panel], ...
               'Units', 'normalized');

% Initialize global variables for filtering, thresholding
  % Load files for analysis and plot in orig_axis ------------------------
  openfiles = UIGetFiles('*.merec', 'Please select file(s):', pwd);
  while length(openfiles) == 0
    openfiles = UIGetFiles('*.merec', 'Please select file(s):', pwd);
  end

% Initialize the 'browser' structure which contains pertinent info about the
% information being viewed by the user
  browser = struct('openfiles', openfiles);
          % 'browser' structure will have elements corresponding to each
          % selected filemod.bandpass = struct;
  for i = 1:1:length(browser)
      browser(i).basename = regexprep(browser(i).openfiles, '.merec', '');
  end     % browser(i).basename has the base name of the input files without
          % the '.merec' extension now
  
  % 'browser' values that are currently set by initial experimental settings
  % Future intent: * collect these data from each file's header
  browser(1).samplefreq = 10000;
  browser(1).thresh_scaling = 6;
  browser(1).thresh = [];   % IMPORTANT: this is an option to be passed to
                            % SNIPPETFILE.M
  browser(1).thresh_init = 'off';
  browser(1).channels = 1; % to-do: make this modular
  browser(1).psniptimes = [];
  browser(1).psnippeaks = [];
  browser(1).sniptimes = [];
  browser(1).snippeaks = [];
  browser(1).mod_data = [];
  browser(1).snipoptions = [];  % IMPORTANT: this is a snippet options structure
                            % for passage via SNIPPETFILE.M
  browser(1).sniprange = [-20 40]   % IMPORTANT: this is an option to be
                                    % passed to snippetfile
  browser(1).outfilename = strcat(browser(1).basename, '.ssnp');
  
  % Initialize the 'mod' structure to contain information about the
  % modified waveforms of the user's current data
  mod.bandpass = struct('ison','false');
  mod.bandpass.range = [1 1];
  mod.bandpass.filter = [];
  mod.bandpass.poles = 2;
  mod.bandpass.samplefreq = 10000;
  mod.det.filter = [];
  mod.polarity = 0;

% ---------------------------------------------------------------------
    
   
   % Load data via MERECMM --------------------------------------------
   for i=1:length(openfiles)
     browser(i).memm = merecmm(browser(i).openfiles, 'contiguous', 1);
   end
          % browser(i).memm contains the merecmm objects (header, channels,
          % etc.)  Using a = browser(*).memm(channel#, :) to load a channel
          % contents into local memory
   browser(1).data = browser(1).memm(1,:); %local, modular copy of data
   
   % Use window info to set envelope parameters -----------------------
   set(orig_axis, 'Units', 'pixels');  % change units to screen pixels
   browser(1).orig_pos = get(orig_axis, 'Position'); % copy position info
   browser(1).envelope_size = browser(1).orig_pos(3); % extract x axis length
   set(orig_axis, 'Units', 'normalized') % re-set axis to normalized units
   
   % Apply enveloping to whole data set -------------------------------
   browser(1).envelope = envelopemem(browser(1).data, ...
             (length(browser(1).data)/(browser(1).envelope_size)));
   browser(1).envelope_xunits = [1/length(browser(1).envelope.min):...
                                 1/length(browser(1).envelope.min):1]...
                                 *(length(browser(1).data)/(browser(1).samplefreq));
   
   % Plot data in orig_axis -------------------------------------------
   browser(1).xpts = 1:1:length(browser(1).envelope.min); %to be removed
   set(orig_axis, 'Xlim', [0 length(browser(1).data)/browser(1).samplefreq]);
   browser(1).whereami = [1 length(browser(1).data)];
   subplot(orig_axis);
   browser(1).envelope_fill = fillmm2(browser(1).envelope.min, ...
                              browser(1).envelope.max, browser(1).envelope_xunits);
   
   % Mirror original data in mod_axis ---------------------------------
   linkaxes([orig_axis mod_axis], 'x');
   browser(1).mod_data = browser(1).data;
   browser(1).mod_envelope = browser(1).envelope;
   browser(1).mod_envelope_xunits = browser(1).envelope_xunits;
   subplot(mod_axis);
   browser(1).mod_fill = fillmm2(browser(1).mod_envelope.min, ...
                   browser(1).mod_envelope.max, browser(1).mod_envelope_xunits);                       
                          
% EMBEDDED FUNCTION DEFINITIONS =======================================
    % Define embedded function 'zoombtn_Callback' for the Zoom button 
    % at the top-right side the orig_panel
    function zoombtn_Callback(source, data)
        prnt = get(source, 'Parent');
        switch prnt
            case orig_panel
                subplot(orig_axis);
                p = pan;
                set(p, 'Enable', 'off'); %Turn PAN mode OFF
                zoom(orig_axis, 'on'); %Turn ZOOM mode ON
        end
    end
    % End of 'zoombtn_Callback defn ---------------------------------
    
    % Define embedded function 'panbtn_Callback' for the Pan button
    % at the top-right side of the orig_panel
    function panbtn_Callback(source, data)
        prnt = get(source, 'Parent');
        switch prnt
            case orig_panel
                zoom(orig_axis, 'off'); % Turn ZOOM mode OFF
                pan(orig_axis,'xon'); % Turn PAN mode to XON (x-only)
                subplot(orig_axis);
                p = pan; % assign a value to the pan mode on this axis
                set(p, 'ActionPostCallback', {@pan_over}); % callback
                            % 'pan_over' called when pan action completes
        end
    end
    % End of 'panbtn_Callback defn ----------------------------------

    % Define embedded function 'envbtn_Callback' for the Envelope
    % button at the top-right of the orig_panel
    function envbtn_Callback(source, data)
       set(orig_axis, 'Units', 'pixels');
       position = get(orig_axis, 'Position'); % get position
       browser(1).envelope_size = position(3); % set envelope
                                        %size to the # pixels on axis
       xlim = (get(orig_axis, 'Xlim')); % assign xlim
                    % as the current time range displayed after any zooming
       if xlim(1)<0
            xlim(1) = 0; % if they zoom past the beginning, set
                         % xlim lower limit at 0
       end
                %scanperpix = floor((browser(1).whereami(2)-browser(1).whereami(1))...
                    %/browser(1).envelope_size);
                    % assign scanperpix as the # of actual data scans
                    % represented within the zoomed region
       xlimscan = int32(xlim*browser(1).samplefreq); % assign xlimscan as the
                    % axis' limits in terms of scan # instead of pixel
                    % number
       if xlimscan(1) == 0
            xlimscan(1) = 1;
       end
       browser(1).whereami = xlimscan;
                    % update 'whereami' to reflect the current scan
                    % position of the axis for use in future
                    % zooming/panning
       newdata = browser(1).data(browser(1).whereami(1)+1:...
                          browser(1).whereami(2));
       browser(1).envelope_xunits = [double(browser(1).whereami(1)):...
                               (double(browser(1).whereami(2)+1-browser(1).whereami(1))/browser(1).envelope_size):...
                                double(browser(1).whereami(2))]/browser(1).samplefreq;
                                   % assign 'newdata' as the data for the scans within the
                                   % zoomed region       
       newmoddata = browser(1).mod_data(browser(1).whereami(1):...
                          browser(1).whereami(2));
       if browser(1).envelope_size > (browser(1).whereami(2)-...
                                      browser(1).whereami(1));
                              % if the envelope size (# pixels in x-axis)
                              % is greater than the number of scans to be
                              % viewed:
          browser(1).envelope_xunits = double([browser(1).whereami(1):...
                        1:browser(1).whereami(2)])/browser(1).samplefreq;
          plot(orig_axis, browser(1).envelope_xunits,...
                         browser(1).data(browser(1).whereami(1):...
                         browser(1).whereami(2)));
          plot(mod_axis, browser(1).envelope_xunits,...
                         browser(1).mod_data(browser(1).whereami(1):...
                         browser(1).whereami(2)));
                              % only plot the raw data without enveloping
          return
       end
       browser(1).envelope = envelopemem(newdata, ...
                    (length(newdata)/(browser(1).envelope_size)));
       browser(1).mod_envelope = envelopemem(newmoddata,...
                    (length(newmoddata)/(browser(1).envelope_size)));
                        % assing browser(1).envelope as the # scans
                        % per upcoming envelope (?? same as xlimscan??)
                % browser(1).xpts = 1:1:length(browser(1).envelope.min);
                        % assign 'xpts' as the x pixel values for
                        % enveloping (req'd for fillmm2 call below)
       cla(orig_axis, 'reset'); % clear the axis, reset limits
       cla(mod_axis, 'reset'); 
       set(orig_axis, 'Xlim', [double(browser(1).whereami(1))/browser(1).samplefreq ...
                               double(browser(1).whereami(2))/browser(1).samplefreq]);
                        % force the axis limits to be the # of pixels of
                        % x-axis
       subplot(orig_axis); % set orig_axis to current axis
       browser(1).envelope_fill = fillmm2(browser(1).envelope.min, ...
                              browser(1).envelope.max, browser(1).envelope_xunits);
                        % use fillmm2 to plot the new envelope
       subplot(mod_axis);
       browser(1).mod_envelope_fill = fillmm2(browser(1).mod_envelope.min, ...
                      browser(1).mod_envelope.max, browser(1).envelope_xunits);
       set(orig_axis, 'Units', 'normalized');
       if strcmp(browser(1).thresh_init,'on')
           mod_threshline = 0;
           thresh_slider_update(source,data);
           if snip_plot
               calc_snips(source,data);
           end
       end
    end
    %end 'envbtn_Callback' defn

    % Define embedded function resetbtn_Callback to re-set the orig_axis
    % to the entire data section
    function resetbtn_Callback(source, data)
        subplot(orig_axis);
        p = pan;
        set(p, 'Enable', 'off'); %Turn PAN mode OFF
        subplot(mod_axis);
        p = pan;
        set(p, 'Enable', 'off'); %Turn PAN mode OFF
        prnt = get(source, 'Parent');
        switch prnt
            case orig_panel
                browser(1).whereami = [1 length(browser(1).data)];
                        % re-set whereami to the full data section
                set(orig_axis, 'Units', 'pixels');
                position = get(orig_axis, 'Position');
                browser(1).envelope_size = position(3);
                        % browser(1).envelope_size is the # pixels
                        % in the x-axis of the plot
                browser(1).envelope = envelopemem(browser(1).data, ...
                    (length(browser(1).data))/(browser(1).envelope_size));
                browser(1).mod_envelope = envelopemem(browser(1).mod_data,...
                    (length(browser(1).data))/(browser(1).envelope_size));
                        % set browser(1).envelope (see envbtn_Callback
                        % above)
                browser(1).envelope_xunits = [1/length(browser(1).envelope.min):...
                                 1/length(browser(1).envelope.min):1]...
                                 *(length(browser(1).data)/(browser(1).samplefreq));
                cla(orig_axis, 'reset');
                cla(mod_axis, 'reset');
                subplot(orig_axis);
                browser(1).envelope_fill = fillmm2(browser(1).envelope.min, ...
                              browser(1).envelope.max, browser(1).envelope_xunits);
                        % use fillmm2 to draw envelope
                subplot(mod_axis);
                browser(1).envelope_fill = fillmm2(browser(1).mod_envelope.min,...
                              browser(1).mod_envelope.max, browser(1).envelope_xunits);
                set(orig_axis, 'Xlim', double(browser(1).whereami)/browser(1).samplefreq);
                set(orig_axis, 'Units', 'normalized');
            end % of SWITCH
         if strcmp(browser(1).thresh_init,'on')
             mod_threshline = 0;    
             thresh_slider_update(source,data);
             if snip_plot
                 calc_snips(source,data);
             end
         end
    end
    % End resetbtn_Callback defn -------------------------------------

    % Define embedded function 'pan_over' to re-draw the intended section
    % using enveloping after a PAN action is competed by user
    function pan_over(obj,env)
        % see envbtn_Callback above for much of the following section:
        set(orig_axis, 'Units', 'pixels')
        position = get(orig_axis, 'Position');
        browser(1).envelope_size = position(3);
        xlim = get(orig_axis, 'Xlim');
        set(orig_axis, 'Units', 'normalized');
        %scanperpix = floor((browser(1).whereami(2)-browser(1).whereami(1))...
        %           /browser(1).envelope_size);
        xlimscan = int32(xlim*browser(1).samplefreq);
        browser(1).whereami = xlimscan;
            if browser(1).envelope_size > (xlimscan(2)-xlimscan(1))
                browser(1).envelope_xunits = double([browser(1).whereami(1):...
                     1:browser(1).whereami(2)])/browser(1).samplefreq;
                plot(orig_axis, browser(1).envelope_xunits,...
                    browser(1).data(browser(1).whereami(1):...
                    browser(1).whereami(2)));
                plot(mod_axis, browser(1).envelope_xunits,...
                    browser(1).mod_data(browser(1).whereami(1):...
                    browser(1).whereami(2)));
    
                % only plot the raw data without enveloping
                    % THIS CURRENTLY CAUSES ERRORS FOR AN UNKNOWN REASON
                    % SEEMS TO RELATE TO USE OF 'PLOT' with
                    % 'ActionPostCallback'
            else 
                browser(1).envelope_xunits = [double(browser(1).whereami(1)):...
                  (double(browser(1).whereami(2)+1-browser(1).whereami(1))/browser(1).envelope_size):...
                   double(browser(1).whereami(2))]/browser(1).samplefreq; 
                newdata = browser(1).data(browser(1).whereami(1):...
                          browser(1).whereami(2));
                newmoddata = browser(1).mod_data(browser(1).whereami(1):...
                          browser(1).whereami(2));
                browser(1).envelope = envelopemem(newdata, ...
                    (length(newdata)/(browser(1).envelope_size)));
                browser(1).mod_envelope = envelopemem(newmoddata,...
                    (length(newmoddata)/(browser(1).envelope_size)));
                %browser(1).xpts = 1:1:length(browser(1).envelope.min);
                cla(orig_axis, 'reset');
                cla(mod_axis, 'reset');
                subplot(orig_axis);
                set(orig_axis, 'Xlim', double(browser(1).whereami)/browser(1).samplefreq);
                browser(1).envelope_fill = fillmm2(browser(1).envelope.min, ...
                          browser(1).envelope.max, browser(1).envelope_xunits);
                subplot(mod_axis);
                browser(1).mod_envelope_fill = fillmm2(browser(1).mod_envelope.min, ...
                          browser(1).mod_envelope.max, browser(1).envelope_xunits);
            end
    end
    % End of pan_over callback defn------------------------
    
    % define embedded function bandpass_apply_Callback to apply filter
    % settings to loaded data
    function bandpass_apply_Callback(source, data)
                if get(bandpass_check,'Value') == 1
                    if str2num(get(bandpass_low,'String'))...
                               &str2num(get(bandpass_high, 'String'))
                        mod.bandpass.filter = [];
                        mod.bandpass.poles = 2;          % I have mis-associated 'poles' with order...now fixed @ 2 without changing the variable name
                        mod.bandpass.range = [str2num(get(bandpass_low, 'String')) ...
                                              str2num(get(bandpass_high, 'String'))];
                        [b, a] = cheby1(mod.bandpass.poles, 0.5,...
                                 mod.bandpass.range/(mod.bandpass.samplefreq/2));
                        mod.bandpass.filter(1,:) = b;
                        mod.bandpass.filter(2,:) = a;
                        browser(1).mod_data = filter(mod.bandpass.filter(1,:),...
                                              mod.bandpass.filter(2,:),browser(1).data);
                        browser(1).mod_data(1:1000)=0; % eliminate transient @ from pt 0 to 1
                        envbtn_Callback(source, data);                  
                    end
                else
                    browser(1).mod_data = browser(1).data;
                end % end of "bandpass filtering" check
    end
    %end of bandpass_apply_Callback defn----------------
    
    %embedded function 'thresh_slider_move'
    function thresh_slider_move(source,data)
       browser(1).thresh_scaling = get(thresholding_slider, 'Value'); % update scaling value
       thresh_slider_label_def = sprintf('The current threshold scaling value is: %0.2f',...
                                 get(thresholding_slider, 'Value'));
       set(thresholding_slider_label, 'String', thresh_slider_label_def); %update label
    end
       
    % embedded function thresh_slider_update
    function thresh_slider_update(source,data)
       browser(1).thresh_scaling = get(thresholding_slider, 'Value'); % update scaling value
       thresh_slider_label_def = sprintf('The current threshold scaling value is: %0.2f',...
                                 get(thresholding_slider, 'Value'));
       set(thresholding_slider_label, 'String', thresh_slider_label_def); %update label
       if thresh_mean_noise == 0
           set(thresh_mean_noise_label, 'String', thresh_mean_noise_def)
       else
           tmn = sprintf('Mean noise: %0.5f', browser(1).thresh_mean_noise);
           set(thresh_mean_noise_label, 'String', tmn);
           tl = sprintf('Threshold: %0.5f', ...
                         browser(1).thresh_scaling*browser(1).thresh_mean_noise);
           set(thresh_value_label, 'String', tl);
           if strcmp(browser(1).thresh_init,'on')
               subplot(mod_axis);
               hold on;
               if mod_threshline
                   set(mod_threshline,'Visible', 'off');
               end
               xlims = get(mod_axis, 'Xlim');
               set(mod_axis, 'Units', 'pixels');
               pix = get(mod_axis, 'Position');
               pix = pix(3);
               set(mod_axis, 'Units', 'normalized');
               mod_threshline = line([xlims(1) xlims(2)],...
                           [browser(1).thresh_scaling*browser(1).thresh_mean_noise ...
                            browser(1).thresh_scaling*browser(1).thresh_mean_noise],...
                            'Color','red');
               [browser(1).psniptimes,browser(1).psnippeaks] = ...
                            pseudosnip(browser(1).mod_data,browser(1).thresh,...
                            browser(1).samplefreq,200);
               subplot(mod_axis);
               psnip_plot = scatter(browser(1).psniptimes,browser(1).psnippeaks,8,'g');
               if browser(1).psniptimes
                   tsnpnum = sprintf('Pseudo-\nsnips: %d', length(browser(1).psniptimes));
                   set(threshold_pseudosnip_num, 'String', tsnpnum);
               end
               hold off;
           end
       end
    end

    % embedded function find_mna to get mean noise level
    function find_mna(source,data)
        if ~browser(1).mod_data
            return
        else
            mn = median(browser(1).mod_data,2);
            wf = browser(1).mod_data - repmat(mn, 1, size(browser(1).mod_data,2));
            mna = mean(abs(wf),2);
            browser(1).thresh_mean_noise = mna;
            thresh_mean_noise = 1;
        end
        browser(1).thresh = [-browser(1).thresh_scaling*mna; ...
                              browser(1).thresh_scaling*mna];
        browser(1).thresh_init = 'on';
        thresh_slider_update(source,data);
    end
    
    % define function pseudosnip to give an idea of which events will be
    % detected (and how many of them for a given threshold choice)
    function [tsnips,psnips] = pseudosnip(waveform, thresh, fs, binsize)
    % tsnip = array of time values for pseudosnippets
    % psnip = array of peak values (pos and neg) for pseudosnippets
    % waveform = a 1-d array of values
    % thresh = a scalar describing the amplitude threshold (pos and neg)
    % binsize = number of samples per bin
        bins = length(waveform)/binsize; % the number of in this waveform bins
        bintimes = [(binsize/2):(binsize):length(waveform)+1]/fs;
        psnips = [];
        tsnips = [];
        xhold = 1;
        for i=1:bins
            range = [(i-1)*binsize+1 i*(binsize)];
            imax = max(waveform(range(1):range(2)));
            imin = min(waveform(range(1):range(2)));
            if (imax>thresh)|(imin>thresh)
                tsnips(xhold) = bintimes(i);
                if abs(imin)>imax
                    psnips(xhold)=imin;
                else
                    psnips(xhold)=imax;
                end
                xhold=xhold+1;
            end
            
        end
        %sniptimes = [1/
    end
    
    % define calc_snips function to use snippetfile to return the actual
    % output of snippetfile for the given settings chosen by user
    function calc_snips(source, data)
        browser(1).snipoptions.tofile = 0;
        browser(1).snipoptions.condfilta = mod.bandpass.filter(2,:);
        browser(1).snipoptions.condfiltb = mod.bandpass.filter(1,:);
        browser(1).snipoptions.detfilt = mod.det.filter;
        if get(polarity_both, 'value') == 1
            browser(1).snipoptions.polarity = 0;
        elseif get(polarity_pos, 'value') == 1
            browser(1).snipoptions.polarity = 1;
        elseif get(polarity_neg, 'value') == 1
            browser(1).snipoptions.polarity = -1;
        else
            browser(1).snipoptions.polarity = 0;
        end
        if str2num(get(thresh_close_box, 'String'))
            browser(1).snipoptions.close = str2num(get(thresh_close_box, ...
                                                    'String'));
        else
            browser(1).snipoptions.close = 30;
        end
        tdepth = str2num(get(thresh_troughdepth_box, 'String'));
        if tdepth & tdepth > 0 & tdepth < 1
                browser(1).snipoptions.troughdepth = tdepth;
        else
                browser(1).snipoptions.troughdepth = 0.5;
        end
        if str2num(get(thresh_peaktrough_box, 'String'))
            browser(1).snipoptions.peaktrough = str2num(...
                                   get(thresh_peaktrough_box, 'String'));
        else
            browser(1).snipoptions.peaktrough = 30;
        end
        rebS = str2num(get(thresh_rebnd_S_box, 'String'));
        if rebS & rebS > 0 & rebS < 2000
            browser(1).snipoptions.reboundS = rebS;
        else
            browser(1).snipoptions.reboundS = 100;
        end
        rebC = str2num(get(thresh_rebnd_C_box, 'String'));
        if rebC & rebC > 0 & rebC < 10000
            browser(1).snipoptions.reboundC = rebC;
        else
            browser(1).snipoptions.reboundC = 4000;
        end
        browser(1).snipoptions.rebound0 = 0;
        browser(1).snipoptions.use_alt_header_file = 1;
        header = readheader(browser(1).openfiles);
        browser(1).digital_thresh = uencode(double(browser(1).thresh), 12,...
                                            header.voltageMax,'signed');  % the 0.5 value needs
                                                              % to map to
                                                              % the header
                                                              % file!!!
        [browser(1).sniptimes, browser(1).snippeaks] = ...
            snippetfile(browser(1).openfiles,...
                        'all',...  %scanrange flag
                        1,... %channel number... in future should be modifiable
                        browser(1).digital_thresh, ... %important: threshold
                        browser(1).sniprange, ... %important: sniprange
                        browser(1).snipoptions);
        browser(1).sniptimes = browser(1).sniptimes{1}/10000;
        snippeakcopy = cell2mat(browser(1).snippeaks);
        snpk = int16(max(snippeakcopy(:,:)));
        browser(1).snippeaks = udecode(snpk,12,header.voltageMax);
        subplot(mod_axis);
        hold on;
        if psnip_plot
            set(psnip_plot, 'Visible', 'off');
            snip_plot = scatter(browser(1).sniptimes(5:length(browser(1).sniptimes)),...
                 browser(1).snippeaks(5:length(browser(1).snippeaks)),8,'c');
        end
        thresh_snip_value = sprintf('The number of\nreal snips is: %d', ...
                                     length(browser(1).sniptimes));
        set(thresh_snip_label, 'String', thresh_snip_value);
    end
    
    % define write_ops function to write a text file of the current options
    % for this file chosen by the user for later use
    function write_ops(source, data)
        snops_file = strcat(browser(1).basename, '.snops');
        time_stamp = datestr(now);
        writefile = fopen(snops_file, 'wt');
        %write stamp and version info:
        fprintf(writefile, strcat(time_stamp, '\n'));
        fprintf(writefile, 'Current merec_gui_thresholding version: %0.1f\n', version_id);
        %write filename
        fprintf(writefile, strcat(browser(1).openfiles, '\n'));
        %write channel info (to-do: needs to be modifiable in future)
        fprintf(writefile, 'channel#\t%d\n', 1);
        %write bandpass data
        fprintf(writefile, 'Bandpass Applied:\t%d\n', get(bandpass_check,'Value'));
        fprintf(writefile, 'Number of Poles:\t%d\n', mod.bandpass.poles);
        fprintf(writefile, 'Lower Limit (Hz):\t%d\n', mod.bandpass.range(1));
        fprintf(writefile, 'Upper Limit (Hz):\t%d\n\n', mod.bandpass.range(2));
        %write snipoptions data
        effs = [];
        for i = 1:length(browser(1).snipoptions.condfilta)
            effs = strcat(effs,'%f\t');
        end
        writefilters = strcat('snipoptions.condfilta\t', effs, '\n');
        fprintf(writefile, writefilters, browser(1).snipoptions.condfilta);
        writefilters = strcat('snipoptions.condfiltb\t', effs, '\n');
        fprintf(writefile, writefilters, browser(1).snipoptions.condfiltb);
        fprintf(writefile, 'snipoptions.detfilt\t%f\n',browser(1).snipoptions.detfilt);
        fprintf(writefile, 'snipoptions.polarity\t%d\n',browser(1).snipoptions.polarity);
        fprintf(writefile, 'snipoptions.close\t%d\n',browser(1).snipoptions.close);
        fprintf(writefile, 'snipoptions.troughdepth\t%d\n',browser(1).snipoptions.troughdepth);
        fprintf(writefile, 'snipoptions.peaktrough\t%d\n',browser(1).snipoptions.peaktrough);
        fprintf(writefile, 'snipoptions.reboundS\t%d\n',browser(1).snipoptions.reboundS);
        fprintf(writefile, 'snipoptions.reboundC\t%d\n',browser(1).snipoptions.reboundC);
        fprintf(writefile, 'snipoptions.rebound0\t%d\n',browser(1).snipoptions.rebound0);
        fprintf(writefile, 'snipoptions.use_alt_header_file\t%d\n', browser(1).snipoptions.use_alt_header_file);
        % write snippetfile arguments
        fprintf(writefile, 'snippetfile:infilename\t%s\n',browser(1).openfiles);
        fprintf(writefile, 'snippetfile:scanrange\t%d\t%d\n', 0, length(browser(1).data));
        nchan = [];
        for i = 1:length(browser(1).channels)
            nchan = strcat(nchan, '%d\t');
        end
        chanwrite = strcat('snippetfile:channels\t',nchan,'\n');
        fprintf(writefile, chanwrite, browser(1).channels);
        nthresh = [];
        for i = 1:size(browser(1).digital_thresh,2)
            nthresh = strcat('%f\t');
        end
        threshwrite = strcat('snippetfile:thresh\t',nthresh,'\n');
        fprintf(writefile, threshwrite, browser(1).digital_thresh(1,:));
        fprintf(writefile, threshwrite, browser(1).digital_thresh(2,:));
        fprintf(writefile, 'snippetfile:varargin:sniprange\t%d\t%d\n', browser(1).sniprange);
        fprintf(writefile, 'snippetfile:varargin:outfilename\t%s\n', strcat(browser(1).basename,'.ssnp'));
        fclose(writefile);
    end

    % define snip_it function to take current values and pass them to
    % snippetfiled to create .ssnp files from the data
    function snip_it(source, data)
        browser(1).snipoptions.tofile = 1;
        browser(1).snipoptions.outfilename = browser(1).outfilename;
        snippetfile(browser(1).openfiles,...
                    'all',...
                    browser(1).channels,...
                    browser(1).digital_thresh,...
                    browser(1).sniprange,...
                    browser(1).snipoptions,...
                    browser(1).outfilename);
    end

    % embedded function generate_vlv_file takes the current data and writes
    % a valve file for that input merec file using 'merec2vlv'
    function gen_vlv_file(source, data)
        outfile = strcat(browser(1).basename, '.vlv');
        merecfiles = browser(1).openfiles;
        timestim_ops.nocheck = 0;
        timestim_ops.feedback = 0;
        timestim_ops.timestim_calc_vs_intended_tolerance = 1.0;
        timestim_ops.wndL = 0.2;
        timestim_ops.wndR = 0.2;
        merec2vlv(outfile, merecfiles, timestim_ops);
    end

    function snip_all_in_dir(source, data)
       a = 0; 
    end
    % END of EMBEDDED FUNCTION DEFNS ===============================================    
    
    
end % end of merec_gui_thresholding
