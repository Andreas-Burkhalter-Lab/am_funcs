% Runs all steps in Ryan_Mouse_Stimulus_Protocol prior to first Flip
% clear all;clc

SPEED_MULTIPLIER = 10;
MAX_FRAMES = 500;                                            %Multiply run time by monitor refresh rate to select select the appropriate number
Angle = 270;                                                   %Degrees
modulateColor = [255 255 255];                                %RGB Change to grating color
phase = 0;                                                    %Degrees
freq = .005;                                                   %Spatial Frequency (Hz)
contrast = 500;                                               %Contrast
sc = 250.0;                                                   %Spatial Constant of Gaussian hull function - the sigma value
nonSymmetric = 0;                                             % 0 makes ellipsoid, 1 makes a symmetric circle
aspectratio = 0;                                            % Ignored if nonSymmetric = 1
disableNorm = 0;                                              % If set to a value of 1, the special multiplicative normalization term normf = 1/(sqrt(2*pi) * sc) will not be applied to the computed gabor.

backgroundColorOffset = [0 0 0 0 ];                          %
radius = inf;                                                %  These are optional arguments
contrastPreMultiplicator = 1;                              %


myScreen = max(Screen('Screens'));
[ windowPtr , rect ] = Screen(myScreen,'OpenWindow');

width = rect(3);
height = rect(4);

[gaborid, gaborrect] = CreateProceduralGabor(windowPtr, width, height, nonSymmetric, backgroundColorOffset, disableNorm, contrastPreMultiplicator);
 
dstRect = [];                                  %No offset
% dstRect = OffsetRect(gratingrect,xc, yc);     %Offset the grating from the center at the start

Screen('DrawTexture', windowPtr, gaborid, [], [], Angle, [], [], [modulateColor], [], kPsychDontDoRotation, [180-phase, freq, sc, contrast, aspectratio, 0, 0, 0]);
% 
% % Perform initial flip to gray background and sync us to the retrace:
% vbl = Screen('Flip', windowPtr);
% count = 0;
% 
% while count < MAX_FRAMES
%    count = count + 1;
%     
%    % Set velocity of movement:
%    shift = count*SPEED_MULTIPLIER;
%     
%    Screen('DrawTexture', windowPtr, gaborid, [], [dstRect], Angle, [], [], [modulateColor], [], kPsychDontDoRotation, [180-phase+shift, freq, sc, contrast, aspectratio, 0, 0, 0]);
%     
%    % Go at normal refresh rate for good looking gabors:
%    Screen('Flip', windowPtr);
%    
%      
% end
% 
% Screen('CloseAll')