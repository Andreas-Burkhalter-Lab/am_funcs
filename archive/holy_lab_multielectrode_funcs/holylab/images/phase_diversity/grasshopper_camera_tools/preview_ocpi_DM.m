function preview_ocpi_DM

  num_triggers = inf;

  vid1 = videoinput('dcam',1);
  triggerconfig(vid1, 'manual')
  set(vid1,'TriggerRepeat',num_triggers);
  set(vid1, 'FramesPerTrigger',1)
  src1 = getselectedsource(vid1);
  src1 = set_camera_features(src1);

  vid2 = videoinput('dcam',2);
  triggerconfig(vid2, 'manual')
  set(vid2,'TriggerRepeat',num_triggers);
  set(vid2, 'FramesPerTrigger',1)
  src2 = getselectedsource(vid2);
  src2 = set_camera_features(src2);

  vidRes = get(vid1, 'VideoResolution');
  nBands = get(vid1, 'NumberOfBands');

  im1 = rand(vidRes(2), vidRes(1));
  im2 = im1;

  hFig = figure('Position',[20 50 1220 850],...
    'Toolbar','none',...
    'Menubar', 'none',...
    'NumberTitle','Off',...
    'Name','Dual Camera Preview');

  h(1) = axes('Position',[0.007 0.23 .485 1]);
  him = image( im1 );
  set(him(1),'CDataMapping','scaled')
  set(him(1),'EraseMode','none');
  set(h(1), 'CLimMode','manual');
  set(h(1),'Visible','off');
  set(h(1),'DataAspectRatio',[1 1 1]);
  set(h(1),'DataAspectRatioMode','manual')

  h(2) = axes('Position',[0.505 0.23 .485 1]);
  him(2) = image( im2 );
  set(him(2),'CDataMapping','scaled')
  set(him(2),'EraseMode','none');
  set(h(2), 'CLimMode','manual');
  set(h(2),'Visible','off');
  set(h(2),'DataAspectRatio',[1 1 1]);
  set(h(2),'DataAspectRatioMode','manual')

  check_box_gray = uicontrol('Style', 'checkbox',...
    'String', 'gray scale',...
    'Units','normalized',...
    'Position',[0.01 0.2 .06 .02]);
  set(check_box_gray,'Callback',{@gray_color, h});

  function gray_color(src, eventData, h)
    if (get(check_box_gray,'Value') == get(check_box_gray,'Max'))
      colormap(gray)
    else
      colormap(jet)
    end
  end

  check_box_autoscale = uicontrol('Style', 'checkbox',...
    'String', 'auto scale',...
    'Units','normalized',...
    'Position',[0.01 0.25 .06 .02]);

  display_scale_min = uicontrol('Style','slider',...
    'Min',0,'Max',66000,...
    'Value',0,...
    'SliderStep',[1e-3 1e-2],...
    'Units','normalized',...
    'Position',[0.01 0.3 0.12 0.02]);

  display_text_min = uicontrol('Style','text',...
    'String',['Min Intensity Val = ' num2str(get(display_scale_min,'Value'))],...
    'Units','normalized',...
    'Position',[0.01 0.315 0.12 0.02]);

  display_scale_max = uicontrol('Style','slider',...
    'Min',0,'Max',66000,...
    'Value',66000,...
    'SliderStep',[1e-3 1e-2],...
    'Units','normalized',...
    'Position',[0.01 0.35 0.12 0.02]);

  display_text_max = uicontrol('Style','text',...
    'String',['Max Intensity Val = ' num2str(get(display_scale_max,'Value'))],...
    'Units','normalized',...
    'Position',[0.01 0.365 0.12 0.02]);

  full_range = [1 1 vidRes];
  current_range = full_range;

  set(him,'ButtonDownFcn',@po_set_range);

  function po_set_range(src,eventData)
    selType = get(hFig,'SelectionType');
    switch selType
      case 'normal'
        current_range = round(getrect);
        is_too_small = current_range < 1;
        current_range(is_too_small) = 1;
        shrink_right = current_range(1:2) + current_range(3:4) - vidRes;
        shrink_right(shrink_right < 0) = 0;
        current_range(3:4) = current_range(3:4) - shrink_right;
      case 'alt'
        current_range = full_range;
    end
  end

  uicontrol('String', 'Close',...
    'Callback', @po_quit_previewing,...
    'Units','normalized',...
    'Position',[0.8 0.025 .10 .05]);

  function po_quit_previewing(src,eventData)
    keep_previewing = false;
  end
 
  DM_initialize_button = uicontrol('String', 'DM initialize',...
    'Callback', @DM_initialize,...
    'Units','normalized',...
    'Position',[0.4 0.025 .10 .05]);
  
  set(DM_initialize_button,'UserData','DM_uninitialized');
  
  %------------------------------------------------

  function DM_initialize(src, eventData)
    
    set(DM_initialize_button,'Enable','off');
    set(DM_initialize_button,'UserData','DM_initialized');
      
    % DM related stuff
    addpath('d:\dmirror'); % the DM initialize seems to work only when in that folder
    dm_initialize % will initialize the dm
    pause(3) % give enough time for initialization
    flat = dm_flat;
    dm_apply(flat); % this will send in the first signal. only after the dm_flat can we continue to send other commands
    %pause(2)

    % h(3) is the DM shape
    h(3) = axes('Position',[0.25 0.1 0.3 0.3]);
    him(3) = image(rand(32,32));
    set(him(3),'CDataMapping','scaled')
    set(him(3),'EraseMode','none');
    set(h(3), 'CLimMode','manual');
    set(h(3),'Visible','off');
    set(h(3),'DataAspectRatio',[1 1 1]);
    set(h(3),'DataAspectRatioMode','manual')
    colorbar

    uicontrol('Style','text', 'String', 'DM shape',...
      'Units','normalized',...
      'Position',[0.32 0.41 .10 .02]);
    
 %--------------------------------------------
 % Defocus adjustment

    dm_defocus = uicontrol('Style','slider',...
      'Min',-0.4,'Max',0.4,...
      'Value',0,...
      'SliderStep',[1e-2 1e-1],...
      'Units','normalized',...
      'Position',[0.59 0.35 0.08 0.02]);

    dm_defocus_text = uicontrol('Style','text',...
      'String',['Defocus = ' num2str(get(dm_defocus,'Value'))],...
      'Units','normalized',...
      'Position',[0.57 0.37 0.12 0.02]);

    set(dm_defocus,'Callback',{@dm_defocus_callback, dm_defocus_text, dm_defocus})

    function dm_defocus_callback(src, eventData, dm_defocus_text, dm_defocus)
      set(dm_defocus_text,'String',['Defocus = ' num2str(get(dm_defocus,'Value'))])
      change_dm_shape(src, eventData);
    end
    
    %---------------------------------------------------------------
    % Tissue Angle Adjustment
    
    tissue_angle = uicontrol('Style','slider',...
      'Min',0,'Max',90,...
      'Value',60,...
      'SliderStep',[1e-2 1e-1],...
      'Units','normalized',...
      'Position',[0.59 0.25 0.08 0.02]);

    tissue_angle_text = uicontrol('Style','text',...
      'String',['tissue angle = ' num2str(get(tissue_angle,'Value'))],...
      'Units','normalized',...
      'Position',[0.57 0.27 0.12 0.02]);

    set(tissue_angle,'Callback',{@tissue_angle_callback, tissue_angle_text, tissue_angle})

    function tissue_angle_callback(src, eventData, tissue_angle_text, tissue_angle)
      set(tissue_angle_text,'String',['tissue angle = ' num2str(get(tissue_angle,'Value'))])
      change_dm_shape(src, eventData);
    end
    
    %-----------------------------------------------------------------
    % Refractive Index adjustment
    
    ref_indx = uicontrol('Style','slider',...
      'Min',1.3,'Max',1.5,...
      'Value',1.33,...
      'SliderStep',[1e-2 1e-1],...
      'Units','normalized',...
      'Position',[0.59 0.15 0.08 0.02]);

    ref_indx_text = uicontrol('Style','text',...
      'String',['Refractive Index = ' num2str(get(ref_indx,'Value'))],...
      'Units','normalized',...
      'Position',[0.57 0.17 0.12 0.02]);

    set(ref_indx,'Callback',{@ref_indx_callback, ref_indx_text, ref_indx})

    function ref_indx_callback(src, eventData, ref_indx_text, ref_indx)
      set(ref_indx_text,'String',['Refractive Index = ' num2str(get(ref_indx,'Value'))])
      change_dm_shape(src, eventData);
    end
    
    function change_dm_shape(src, eventData)
      dm_defocus_value = get(dm_defocus,'Value');
      tissue_angle_value = get(tissue_angle,'Value');
      ref_indx_value = get(ref_indx, 'Value');
      
      % [dm_shape, vol] = ocpi_deformation(dm_defocus_value,...
      % tissue_angle_value, ref_indx);
      
      vol = flat + randn(1,52)/50;
      dm_shape = rand(32,32);
      
      set(him(3),'CData',dm_shape);
      dm_apply(vol);
      
    end

  end

  %------------------------------------------------------------------
  %------------------------------------------------------------------
  keep_previewing = true;
  start(vid1);
  start(vid2);

  while keep_previewing

    trigger(vid1);
    trigger(vid2);

    show_range = current_range + [0 0 current_range(1:2)];

    axes(h(1));
    set(him(1),'CData',im1);
    set(h(1),'XLim',show_range([1 3]),'YLim',show_range([2 4]));

    axes(h(2));
    set(him(2),'CData',im2);
    set(h(2),'XLim',show_range([1 3]),'YLim',show_range([2 4]));

    if (get(check_box_autoscale,'Value') == get(check_box_autoscale,'Max'))
      clims(1) = min([min(im1(:)),min(im2(:))]); %minimum of both the cameras
      clims(2) = max([max(im1(:)), max(im2(:))]); %maximum of both the cameras
    else
      clims(1) = round(get(display_scale_min,'Value'));
      clims(2) = round(get(display_scale_max,'Value'));
      if clims(1) == 66000 % to prevent out of range
        clims(1) = 66000 - 1;
      end
      if clims(1) >= clims(2) % to prevent max being less than min
        clims(2) = clims(1) + 1;
        set(display_scale_max, 'Value',clims(2))
      end
    end

    set(h(1), 'CLim',clims);
    set(h(2), 'CLim',clims);
    set(display_text_min,'String',['Min Intensity Val = ' num2str(clims(1))]);
    set(display_text_max,'String',['Max Intensity Val = ' num2str(clims(2))]);

    im1 = getdata(vid1,1);
    im2 = getdata(vid2,1);

  end

  if isequal(get(DM_initialize_button,'UserData'),'DM_initialized')% DM has been initialized
    % uninitialize dm here
    dm_apply(flat);
    dm_cleanup;
  end
  
  delete(vid1); delete(vid2);
  close(hFig);

end


  