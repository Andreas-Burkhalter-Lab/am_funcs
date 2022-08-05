function slidwincmenu(action,parent)
% SLIDWINCMENU: set up context menu for calling sliderwindow in current axis 
% Syntax:
%   slidwincmenu(hax)
% where hax is the list of axis handles you want passed to sliderwindow,
% in the proper order.
  if ishandle(action)
    hax = action;
    if (nargin > 1)
      hobj = parent;
    else
      hobj = gca;
    end
    % Define the context menu and associate it with current axis
    cmenu = get(hobj,'UIContextMenu');  % if already exists, just add on
    if isempty(cmenu)
      cmenu = uicontextmenu;
      set(hobj,'UIContextMenu', cmenu);
    end
    % Define callback for context menu items
    item = uimenu(cmenu, 'Label', 'SliderWindow', ...
                  'Callback', 'slidwincmenu make', ...
                  'UserData', hax);
  
  elseif strcmp(action,'make')
    % Take care of callbacks
    hax = get(gcbo,'UserData');
    hswfig = sliderwindow(hax,struct('restoreonclose',1));
    hswfigstr = num2str(hswfig);
    % Set it so that deletion of any of the axes under control of the
    % slider window causes the slider window figure to be closed
    set(hax,'DeleteFcn',['if ishandle(',hswfigstr,'), close(', ...
                         num2str(hswfig),'), end']);
  end
