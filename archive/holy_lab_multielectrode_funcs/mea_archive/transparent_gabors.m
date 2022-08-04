% % % transparent_gabors
% last edited 4/6/15
% draw a 2 transparent gabors with individual orientations, sf, and tf, with one as an
% annulus around the other

rf_center = [305 284];

%% Parameters
% For all parameters with more than one argument, the first is for rf
% center, the next (if applicable) for center-annulus spacing, and third
% for annulus.
nframes= 100;   %Use 1500 for ~25s stim; Multiply run time by monitor refresh rate to select select the appropriate number
ngabors = 3;    %% must match number of columens of arguments for following parameters
speed = [24.7219 0 5];
Angle = [170 0 270];            % degrees
modulateColor = [[0;255;0;0], [-0;25;-100;255], [100;0;0;-255]];  %RGB Change to grating color
backgroundColorOffset = [[0;0;0;1], [0;0;0;0], [0;0;0;0]];   
phase = [0 0 90];   %Degrees
freq = [0.0059 0.0001 0.01];     %Spatial Frequency - adjusts speed so that temporal frequency does not change
contrast = [inf inf inf]; % 1000;     %Contrast; set to inf for no gradient
sc = [50 50 50];         %Spatial Constant of Gaussian hull function - the sigma value
nonSymmetric = [0 0 0];    % 0 makes ellipsoid, 1 makes a symmetric circle
aspectratio = [100 2 2];         % Ignored if nonSymmetric = 1
disableNorm = [0 0 0];     % If set to a value of 1, the special multiplicative normalization term normf = 1/(sqrt(2*pi) * sc) will not be applied to the computed gabor.
contrastPremultiplicator = [1 1 1];  
rect_width = [2000 2500 3000]; %% width of stim rectangles
rect_height = [2000 2500 3000]; %% height of stim rectangles
background_color = [0 0 0 1];     %% only has an effect if window not already open

if isempty(Screen('Windows'))  5 %% if a window is not open, open a window and configure as per ProceduralGarboriumDemo
    myScreen = max(Screen('Screens'));
    [windowPtr , rect] = PsychImaging('OpenWindow',myScreen,background_color);
    PsychImaging('PrepareConfiguration');
    ifi = Screen('GetFlipInterval', windowPtr);       %% get frame duration (inter-flip interval)
    Screen('BlendFunction', windowPtr, GL_ONE, GL_ONE);
else
    windowPtr = Screen('Windows'); 
    windowPtr = windowPtr(1);
    rect = Screen('rect',windowPtr);
end
    
% rf_center =  follow_cursor;  %% HideCursor;   %% call function to select RF center 

gabortex = CreateProceduralGabor(windowPtr,rect_width(1),rect_height(1),nonSymmetric(1),...
    backgroundColorOffset(:,1),disableNorm(1),contrastPremultiplicator(1));
Screen('DrawTexture', windowPtr, gabortex, [], [], [], [], [], [], [], kPsychDontDoRotation,...
    [phase(1), freq(1), sc(1), contrast(1), aspectratio(1), 0, 0, 0]);
texrect = Screen('Rect', gabortex);
inrect = [texrect'; texrect'];

for i = 1:ngabors  %% center gabors on same point
    destinationRect(:,i) = CenterRectOnPoint([1000 0 rect_width(i) rect_height(i)],rf_center(1),rf_center(2));
end

mypars = [freq; sc; contrast; aspectratio; zeros(3,ngabors)]; % create gabor param matrix for all params except phase

% rect_WH = 2*[rect(3),rect(4)]; % 

width = 16000; %width = rect(3);    %% controls elipsoid shape
height = 16000; %height = rect(4);     %% controls elipsoid shape

% Initially sync us to VBL at start of animation loop..... as per ProceduralGarboriumDemo 
vbl = Screen('Flip', windowPtr);
tstart = vbl;
frame = 0;
   
%% Run 
for frame = 1:nframes 
   Screen('DrawTextures', windowPtr, gabortex, [], destinationRect, Angle, [],... % draw all three gabors simultaneously
        [],[modulateColor], [], kPsychDontDoRotation,[frame*speed+180-phase; mypars]);
   Screen('Flip', windowPtr); 
   Screen('DrawingFinished', windowPtr);      %% recommended by ProceduralGarboriumDemo
end

% Screen('CloseAll')


