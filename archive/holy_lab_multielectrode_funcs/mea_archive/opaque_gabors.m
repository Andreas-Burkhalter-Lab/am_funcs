function [] = opaque_gabors()
% purpose (not yet working): draw multiple overlapping gabors that mask each other


rf_center = [500 200];
% rf_center = follow_cursor;

%% Parameters
% For all parameters with more than one argument, the first is for rf center, 
% the next for annulus, and the third (if applicable) for  center-annulus spacing. 
nframes= 100;   %Use 1500 for ~25s stim; Multip ly run time by monitor refresh rate to select select the appropriate number
ngabors = 3;    %% must match number of columens of arguments for following parameters
speed = [5 2 0];
Angle = [170 0 270];            % degrees
modulateColor = [[255;255;255;0], [0;255;255;100], [100;0;0;0]];  %RGB Change to grating color
backgroundColorOffset = [[0;0;0;0], [0;0;0;0], [0;0;0;0]];   
phase = [0 0 90];   %Degrees
freq = [0.0059 0.001 0.001];     %Spatial Frequency - adjusts speed so that temporal frequency does not change
contrast = [1e16 1e16 1e16]; % 1000;     %Contrast; set to inf for no gradient
sc = [50 50 50];         %Spatial Constant of Gaussian hull function - the sigma value
nonSymmetric = [0 0 0];    % 0 makes ellipsoid, 1 makes a symmetric circle
aspectratio = [2 2 2];         % Ignored if nonSymmetric = 1
disableNorm = [0 0 0];     % If set to a value of 1, the special multiplicative normalization term normf = 1/(sqrt(2*pi) * sc) will not be applied to the computed gabor.
contrastPremultiplicator = [1 1 1];  
gabor_width = [3000 4000 3500]; %%
gabor_height = [3000 4000 3500]; %%
background_color = [0 0 0];     %% only has an effect if window not already open

if isempty(Screen('Windows')) 
    myScreen = max(Screen('Screens'));
    [ windowPtr , rect ] = Screen(myScreen,'OpenWindow',background_color);
else
    windowPtr = Screen('Windows'); 
    windowPtr = windowPtr(1);
    rect = Screen('rect',windowPtr);
end

rect_WH = 2*[rect(3),rect(4)]; % rect_WH = [1000 500]; %% width and height of stim rectangle
destinationRect = CenterRectOnPoint([0 0 rect_WH],rf_center(1),rf_center(2));

% Create objects to draw and draw the first frame
for i = 1:2
    [gaborid(i), gaborrect(i,:)] = CreateProceduralGabor(windowPtr, gabor_width(i), gabor_height(i),... 
        nonSymmetric(i), backgroundColorOffset(:,i), disableNorm(i), contrastPremultiplicator(i));
    Screen('DrawTexture', windowPtr, gaborid(i), [], destinationRect, Angle(i), [], [],...
        modulateColor(:,i),[], kPsychDontDoRotation,...
        [180-phase(i), freq(i), sc(i), contrast(i), aspectratio(i), 0, 0, 0]);
end

%% Run
for frame = 1:nframes 
    for i = 1:2
        Screen('DrawTexture', windowPtr, gaborid(i), [], destinationRect, Angle(i), [], [],...
           modulateColor(:,i), [], kPsychDontDoRotation,...
           [180-phase(i)+frame*speed(i),freq(i), sc(i), contrast(i), aspectratio(i), 0, 0, 0]);
    end
    Screen('Flip', windowPtr); 
end

Screen('CloseAll')
