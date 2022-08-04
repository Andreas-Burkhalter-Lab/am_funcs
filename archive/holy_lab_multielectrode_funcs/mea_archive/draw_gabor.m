%% draw_gabor
% last edit 5/4/15
function [] = draw_gabor()
    Screen('Preference', 'SkipSyncTests', 1)
    debugoriginal = Screen('Preference','VisualDebugLevel');
    Screen('Preference', 'VisualDebugLevel', 1);    % set splash screen color to black... comment out for experiments

    speed = 24.7219;
    nframes= 100;        %Use 1500 for ~25s stim; Multiply run time by monitor refresh rate to select select the appropriate number
    Angle = 170;            %Degrees
    modulateColor = [255 255 255 0];  %RGB Change to grating color
    backgroundColorOffset = [0 0 0 0];   
    phase = 0;   %Degrees
    freq = .0059;     %Spatial Frequency - adjusts speed so that temporal frequency does not change
    contrast = 1e38; % % Gabor becomes invisible if this value =>1e(38.6)
    sc = 50.0;         %Spatial Constant of Gaussian hull function - the sigma value
    nonSymmetric = 0;    % 0 makes ellipsoid, 1 makes a symmetric circle
    aspectratio = 1;         % Ignored if nonSymmetric = 1
    disableNorm = 0;     % If set to a value of 1, the special multiplicative normalization term normf = 1/(sqrt(2*pi) * sc) will not be applied to the computed gabor.
    radius = 0.1;         %  These are optional arguments
    contrastPremultiplicator = 1;  
    background_color = [0 0 0];     %% only works if window not already open

    if isempty(Screen('Windows')) 
        myScreen = max(Screen('Screens'));
        [ windowPtr , rect ] = Screen(myScreen,'OpenWindow',background_color);
    else
        windowPtr = Screen('Windows'); 
        windowPtr = windowPtr(1);
        rect = Screen('rect',windowPtr)
    end

%     rf_center = [1000 500];
      rf_center = rect(3:4)/2;
%     rf_center =  follow_cursor;  %% HideCursor;   %% call function to select RF center 

    rect_WH = [1680 1680];  % proportional to Gabor size
%     rect_WH = 2*[rect(3),rect(4)];  %% controls width and height of stim rectangle
    destinationRect = CenterRectOnPoint([0 0 rect_WH],rf_center(1),rf_center(2))
%     destinationRect = []; 
    width = 2000; %width = rect(3);    %% controls ellipsoid shape
    height = 2000; %height = rect(4);     %% controls ellipsoid shape

    [gaborid, gaborrect] = CreateProceduralGabor(windowPtr, width, height, nonSymmetric,...
        backgroundColorOffset, disableNorm, contrastPremultiplicator);
    Screen('DrawTexture', windowPtr, gaborid, [], destinationRect, Angle, [], [],...
        [modulateColor],[], kPsychDontDoRotation,...
        [180-phase, freq, sc, contrast, aspectratio, 0, 0, 0]);
    vbl = Screen('Flip', windowPtr);

    for frame = 1:nframes 
       Screen('DrawTexture', windowPtr, gaborid, [], destinationRect, Angle, [], [],...
           [modulateColor], [], kPsychDontDoRotation,...
           [180-phase+frame*speed,freq, sc, contrast, aspectratio, 0, 0, 0]);
       Screen('Flip', windowPtr); 
    end

%     Screen('CloseAll')
    Screen('Preference', 'SkipSyncTests', 0);
    if exist('debugoriginal','var')
        Screen('Preference', 'VisualDebugLevel', debugoriginal);
    end

end