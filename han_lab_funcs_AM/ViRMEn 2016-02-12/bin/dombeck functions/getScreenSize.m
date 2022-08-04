function scr = getScreenSize()

monitorPositions = get(0, 'MonitorPositions');

%%define virtual screen size based on experimental set up

scr = monitorPositions(1,:);
