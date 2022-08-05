function fmb_moviewatch(mov, fmb, varargin)
% (script/fxn in progress to allow user to play movie and watch/align .fmb
% data

%% handle varargin
if nargin > 2
    if isnumeric(varargin{1})
        orientation = varargin{1};
    else
        orientation = 90;
    end
else
    orientation = 90;
end


%% load data
% movie (if not already loaded)
%mov = yuv2mov('prg002/mov003.yuv', 352,240,'420');

% fmb data
%load('2008_06_04_0123overnight.fmbn', '-mat');

diff_fmb = zeros(4, size(diff(fmb.combined_nonsparse_events(1,:)),2));
diff_fmb(1,:) = diff(fmb.combined_nonsparse_events(1,:));
diff_fmb(2,:) = diff(fmb.combined_nonsparse_events(2,:));
diff_fmb(3,:) = diff(fmb.combined_nonsparse_events(3,:));
diff_fmb(4,:) = diff(fmb.combined_nonsparse_events(4,:));

fmbt = fmb.nonsparse_times;
%% set up axes
if orientation == 0
    mainfig = figure('units', 'pixels', 'position', [100 16 704 520]);
elseif orientation == 90
    mainfig = figure('units', 'pixels', 'position', [100 16 704 704]);
end

temphax = axes('parent', mainfig, 'units', 'normalized', 'position', [0.1 0.1 0.8 0.8]);
temphax = SplitVert([0.5 0.52], [1 0 1], temphax);
movax = temphax(1);
% movax = SplitHoriz([0.37 0.4], [1 0 1], temphax(1));
% panel_axis = ghostaxis(movax(1));
% set(movax(1), 'visible', 'off');
% movax = movax(2);
set(movax, 'xtick', [], 'ytick', [], 'box', 'off', 'color', [0.8 0.8 0.8], 'visible', 'off');
temphax = SplitVert([0.45 0.55], [1 0 1], temphax(2));
tophax = SplitHoriz([0.45 0.55], [1 0 1], temphax(1));
bothax = SplitHoriz([0.45 0.55], [1 0 1], temphax(2));
ax1 = tophax(1);
ax2 = bothax(1);
ax3 = tophax(2);
ax4 = bothax(2);

% put some buttons on there
% figure for control:
controlfig = figure('units', 'pixels', 'position', [100 286 204 240]);
set(get(controlfig,'children'), 'visible', 'off');
remote_control_panel = uipanel(controlfig, 'title', 'Playback Controls', 'units', 'normalized',...
                               'position', [0.02 0.02 0.96 0.96], 'backgroundcolor', [0.8 0.8 0.8]);
play_btn = uicontrol(remote_control_panel, 'Units', 'normalized',...
                          'Position', [0.35 0.8,0.25,0.180], 'String', '>',...
                          'style', 'pushbutton', 'tooltipstring', 'PLAY',...
                          'callback', @playbtn);
setappdata(movax, 'multiplier', 1);
setappdata(movax, 'direction', 1);
pause_btn = uicontrol(remote_control_panel, 'Units', 'normalized',...
                          'Position', [0.35 0.6,0.25,0.180], 'String', '||',...
                          'style', 'pushbutton', 'tooltipstring', 'PAUSE',...
                          'callback', @pausebtn);
ffwd_btn = uicontrol(remote_control_panel, 'Units', 'normalized',...
                          'Position', [0.69 0.8,0.25,0.180], 'String', '>>',...
                          'style', 'pushbutton', 'tooltipstring', 'FAST FWD',...
                          'callback', @ffwdbtn);
rwd_btn = uicontrol(remote_control_panel, 'Units', 'normalized',...
                          'Position', [0.02 0.8,0.25,0.180], 'String', '<<',...
                          'style', 'pushbutton', 'tooltipstring', 'REWIND',...
                          'callback', @rewindbtn);
offset_edit = uicontrol(remote_control_panel, 'Style', 'edit', ...
               'Units', 'normalized', 'Position', [0.52,0.42,0.46,0.16],...
               'String', '-5','callback', @setoffset); % use str2num to get it
offset_label = uicontrol(remote_control_panel, 'Style', 'text', ...
               'Units', 'normalized', 'Position', [0.02,0.42,0.46,0.16],...
               'String', 'Movie offset (frames)');  
setappdata(movax, 'offset', -5);
current_edit = uicontrol(remote_control_panel, 'Style', 'edit', ...
               'Units', 'normalized', 'Position', [0.52,0.24,0.46,0.16],...
               'String', '1', 'callback', @setcurframe); % use str2num to get it
current_label = uicontrol(remote_control_panel, 'Style', 'text', ...
               'Units', 'normalized', 'Position', [0.02,0.24,0.46,0.16],...
               'String', 'Current frame:'); % use str2num to get it
setappdata(movax, 'curframe', 1);
snailbtnpos = uicontrol(remote_control_panel, 'Units', 'normalized',...
                          'Position', [0.63 0.6,0.1,0.180], 'String', '>',...
                          'style', 'pushbutton', 'tooltipstring', 'ONE FRAME FORWARD',...
                          'callback', @snailposbtn);
snailbtnneg = uicontrol(remote_control_panel, 'Units', 'normalized',...
                          'Position', [0.23 0.6,0.1,0.180], 'String', '<',...
                          'style', 'pushbutton', 'tooltipstring', 'ONE FRAME BACK',...
                          'callback', @snailnegbtn);
slowbtnpos = uicontrol(remote_control_panel, 'Units', 'normalized',...
                          'Position', [0.75 0.6,0.15,0.180], 'String', '|>',...
                          'style', 'pushbutton', 'tooltipstring', 'SLOW FORWARD',...
                          'callback', @slowposbtn);
slowbtnneg = uicontrol(remote_control_panel, 'Units', 'normalized',...
                          'Position', [0.06 0.6,0.15,0.180], 'String', '<|',...
                          'style', 'pushbutton', 'tooltipstring', 'SLOW BACK',...
                          'callback', @slownegbtn);
setappdata(movax, 'framerate', 100);

%% Start it all up!

%orient movie for viewer
%fprintf('rotating movie ...');
%if orientation == 90
%   for idx = 1:size(mov,2)
%       for rgbidx = 1:size(mov(idx).cdata, 3)
%           temp(rgbidx) = rot90(mov(idx).cdata(:,:,rgbidx), -1);
%       end
%       mov(idx).cdata = temp;
%   end
%end
%fprintf('done!\n');

% master idx
movie_maxsize = size(mov,2);
% play movie
hold on;
toplot = diff_fmb;
% frametext = text('parent', movax, 'string', ['frame' num2str(1)],...
%     'position', [0.3 0.1]);
idx = 1;
while idx >= 1 && idx <= movie_maxsize
    framerate = getappdata(movax, 'framerate');
    movframe = getappdata(movax, 'curframe');
    multiplier = getappdata(movax, 'multiplier');
    direction = getappdata(movax, 'direction');
    offset = getappdata(movax, 'offset');
    idx = movframe-offset;
    
    if idx+offset <= 0
        movframe = 1;
    else
        movframe = idx+offset;
    end
    
    movie(movax, mov, [1 movframe], framerate, [146 2 0 0]);

    if idx-20 < 1
        padding = zeros(1,200-idx*10);
        time_toplot = [padding fmbt(1, 1:idx*10)];
        ax1_toplot = [padding toplot(1,1:idx*10)];
        ax2_toplot = [padding toplot(2,1:idx*10)];
        ax3_toplot = [padding toplot(3,1:idx*10)];
        ax4_toplot = [padding toplot(4,1:idx*10)];
    else
        time_toplot = fmbt(1,idx*10-200:idx*10);
        ax1_toplot = toplot(1,idx*10-200:idx*10);
        ax2_toplot = toplot(2,idx*10-200:idx*10);
        ax3_toplot = toplot(3,idx*10-200:idx*10);
        ax4_toplot = toplot(4,idx*10-200:idx*10);
    end
    
    plot(time_toplot/1000, ax1_toplot,  'parent', ax1);
    plot(time_toplot/1000, ax2_toplot,  'parent', ax2);
    plot(time_toplot/1000, ax3_toplot,  'parent', ax3);
    plot(time_toplot/1000, ax4_toplot,  'parent', ax4);
    set([ax1 ax2 ax3 ax4], 'xlim', [time_toplot(1)/1000 time_toplot(end)/1000]);
    
    
    if idx+direction*multiplier < 1 || idx+direction*multiplier > movie_maxsize
        pause on;
        pause;
    end
    if movframe ~= getappdata(movax, 'curframe')
        idx = getappdata(movax, 'curframe')-offset;
    end
    idx = idx+(direction)*(multiplier);
    movframe = idx + offset;
    setappdata(movax, 'curframe', movframe);
%    set(frametext, 'string', ['frame' num2str(movframe)]);
    set(current_edit, 'string', num2str(movframe));
    
    %pause(0.01);
end
% play the axes

function playbtn(hObject, eventdata)
   setappdata(movax, 'direction', 1);
   setappdata(movax, 'multiplier', 1);
   setappdata(movax, 'framerate', 100);
end
function pausebtn(hObject, eventdata)
   setappdata(movax, 'multiplier', 0);
end
function setcurframe(hObject, eventdata)
   newval = str2num(get(current_edit, 'string'));
   if ~isempty(newval)
       if newval >= 1 && newval <= movie_maxsize
           setappdata(movax, 'curframe', newval);
       end
   else
       oldval = getappdata(movax, 'curframe');
       set(current_edit, 'string', num2str(oldval));
   end
end
function rewindbtn(hObject, eventdata)
    current = getappdata(movax, 'multiplier');
    if current<=0
        setappdata(movax, 'multiplier', 1);
    else
        setappdata(movax, 'multiplier', current+1);
    end
    setappdata(movax, 'direction', -1);
    setappdata(movax, 'framerate', 100);
end
function ffwdbtn(hObject, eventdata)
    current = getappdata(movax, 'multiplier');
    if current<=0
        setappdata(movax, 'multiplier', 1);
    else
        setappdata(movax, 'multiplier', current+1);
    end
    setappdata(movax, 'direction', 1);
    setappdata(movax, 'framerate', 100);
end
function setoffset(hObject, eventdata)
    current = getappdata(movax, 'offset');
    new = str2num(get(offset_edit, 'string'));
    if ~isempty(new)
        if new ~= current;
            setappdata(movax, 'offset', new);
        end
    end
end
function snailposbtn(hObject, eventdata)
    current = getappdata(movax, 'curframe');
    if current+1 <= movie_maxsize
        setappdata(movax, 'curframe', current+1);
    end
    setappdata(movax, 'multiplier', 0);
end
function snailnegbtn(hObject, eventdata)
    current = getappdata(movax, 'curframe');
    if current+1 <= movie_maxsize
        setappdata(movax, 'curframe', current-1);
    end
    setappdata(movax, 'multiplier', 0);
end
function slowposbtn(hObject, eventdata)
    cur_rate = getappdata(movax, 'framerate');
    cur_dir = getappdata(movax, 'direction');
    if cur_rate >= 100 || cur_rate == 0 || cur_dir == -1
        setappdata(movax, 'framerate', 10);
    else
        setappdata(movax, 'framerate', cur_rate+10);
    end
    setappdata(movax, 'direction', 1);
    setappdata(movax, 'multiplier', 1);
end
function slownegbtn(hObject, eventdata)
    cur_rate = getappdata(movax, 'framerate');
    cur_dir = getappdata(movax, 'direction');
    if cur_rate >= 100 || cur_rate == 0 || cur_dir == 1
        setappdata(movax, 'framerate', 10);
    else
        setappdata(movax, 'framerate', cur_rate+10);
    end
    setappdata(movax, 'direction', -1);
    setappdata(movax, 'multiplier', 1);
end
        

end