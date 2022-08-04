function cursor_position_xy = follow_cursor(varargin)
%FOLLOW_CURSOR: overlays adjustable oriented rectangle on top of cursor
%   [cursorposition] = follow_cursor() outputs the position of the cursor in
%            psychtoolbox global x,y coordinates at termination time. 
% To quit and close the PTB window, use the Escape key, not ctr+c, or else 
%   HideCursor might be left on. 
%%% Last updated 2018/08/28 on thermaltake

width = 50; % bar width in pixels
length = 250; % bar length in pixels
Angle = 30; % orientation of bar in degrees; positive values rotate clockwise

verts = getVerts(Angle, length, width);

% % % bar_dim = [0 0 50 250]; % for fixed-orientation bar
if isempty(Screen('Windows'))
    myScreen = max(Screen('Screens'));
    [ windowPtr , rect ] = Screen(myScreen,'OpenWindow',BlackIndex(myScreen));
else
    windowPtr = Screen('Windows');
    windowPtr = windowPtr(1);
    rect = Screen('Rect', windowPtr);
    Screen('FillRect',windowPtr, BlackIndex(windowPtr), rect); 
end

bar_color = WhiteIndex(windowPtr);
go = true;
HideCursor;
SetMouse(rect(3)/2, rect(4)/2, windowPtr);

%%% key codes from http://www.cambiaresearch.com/articles/15/javascript-char-codes-key-codes
enterCode = 13;%[zeros(1,12), 1 zeros(1,256-13)];
escapeCode = 27; 
leftArrowCode = 37;%[zeros(1,36), 1, zeros(1,256-37)];
upArrowCode = 38;%[zeros(1,37), 1, zeros(1,256-38)]; 
rightArrowCode = 39;%[zeros(1,38), 1, zeros(1,256-39)];
downArrowCode = 40;%[zeros(1,39), 1, zeros(1,256-40)]; 
shiftLeftCode = 160; % really [16 160]
shiftRightCode = 161; % really [16 161]


rotatePerPress = 1; % degrees
lengthChangePerPress = 5; % pixels

% The following lines assumes that the open window is on the non-primary
% screen to the right of the main screen in a multi-display setup. The
% cursor will be centered in the open window.
res = Screen('Resolution',2);
SetMouse(res.width + rect(3)/2, rect(4)/2); 
i = 0; %% initialize counter for vbl checking

fprintf('Press Shift to output current x and y coordinates of bar center on screen.')
while KbCheck; end  % if key is pressed right after function is called, wait for key release
while go
    [xNow yNow] = GetMouse(windowPtr);
    Screen('FillPoly',windowPtr,bar_color,verts+repmat([xNow yNow],size(verts,1),1),1); % oriented bar
% % %     Screen('FillRect',windowPtr,bar_color,CenterRectOnPoint(bar_dim,x,y)); %%% fixed-orientatiom bar 

%%% Use next two lines for the third to check for missed flips.
% % % % % % % % %             i = i+1;
% % % % % % % % %             vbl(i) = Screen('Flip',windowPtr);
    Screen('Flip',windowPtr);

    [junk junk keyCodes] = KbCheck;
    [junk pressed] = find(keyCodes);
    switch max([0 pressed]) % must be scalar
        case shiftLeftCode %|| shiftCodeRight
            break
        case shiftRightCode
            break
        case escapeCode
            SetMouse(res.width/2, res.height/2); % put mouse back in center of command screen
            ShowCursor(windowPtr);
            Screen('Close',windowPtr);
            error('Quit follow_cursor.m')
        case leftArrowCode
            Angle = Angle - rotatePerPress;
            verts = getVerts(Angle, length, width);
        case rightArrowCode
            Angle = Angle + rotatePerPress;
            verts = getVerts(Angle, length, width);
        case downArrowCode
            length = length - lengthChangePerPress;
            verts = getVerts(Angle, length, width);
        case upArrowCode
            length = length + lengthChangePerPress;
            verts = getVerts(Angle, length, width);
    end
end

Screen('FillRect', windowPtr, BlackIndex(windowPtr)) % remove the bar, set screen to black
Screen('Flip', windowPtr)
SetMouse(res.width/2, res.height/2); % put mouse back in center of command screen
ShowCursor(windowPtr);
% Screen(windowPtr,'Close');  %% could also use Screen('MATLABToFront')
cursor_position_xy = [xNow yNow];

end






%% Subfunction for calculating bar vertices based on Angle, length, and width
function [verts] = getVerts(Angle, length, width)
    halfdiag = sqrt((width/2)^2 + (length/2)^2);
    Anglerad = deg2rad(Angle);
    theta = atan(width/length); % corner-center-corner angle in radians
    a1 = -Anglerad+theta; % first corner-center-vertical angle in radians
    a2 = -Anglerad-theta; % second corner-center-vertical angle in radians
    x = halfdiag*sin([a1 a2 a1 a2]).*[1 1 -1 -1]; % vertex x coordinates before centering
    y = halfdiag*cos([a1 a2 a1 a2]).*[1 1 -1 -1]; % vertex y coordinates before centering
    verts = [[x';x(1)], [y';y(1)]];
end