function code = Track_180cm_eh_test1_2Xrew_world_switchRotGain_testing
%modified to remove DAQ and rewards for testing on desk computers
%switches rotation gain in track in new world
% quakeClone1   Code for the ViRMEn experiment quakeClone1. 
%   code = quakeClone1   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.

% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT




% --- INITIALIZATION code: executes before the ViRMEn engine starts.
function vr = initializationCodeFun(vr)
vr.logname=sprintf('Log_%s.dat',datestr(now,30));
vr.fid = fopen(vr.logname,'w');%170512.for logging data
% vr.fid = fopen('virmenLog.dat','w');%170512.for logging data
% vr = initializeDAQ(vr);
vr.timeSolenoid = 140; %in milliseconds
vr.topTarget = 170;%from "linearTrackTwoDirections"
vr.bottomTarget = 10;%new term
% %vr.endZone = 170;  %original
vr.beginZone = 8.5;
vr.end_gain = 2;
vr.track_gain = .4;
vr.RotGainFactor= 1.50;%for changing rot gain in track
vr.RotGainFactorEnd= 0.75;%for changing rot gain in track
vr.currentGoal = 0;
%vr.scaling = [eval(vr.exper.variables.scalingGoal) eval(vr.exper.variables.scalingStart)];
%vr.topTarget = (eval(vr.exper.variables.numCylinders)-1)*eval(vr.exper.variables.cylinderSpacing);
%vr.scalingDecay = eval(vr.exper.variables.scalingDecay);
%vr.currentGoal = 1;
vr.numRewards = 0;%actual earned rewards
vr.water = 0;%keep track of total rewards for water tracking
vr.currentGain = 1;%gain world. 1=default, 2=gainFactor
vr.testRew=0;%track rew num for testing reward amount. see below

% This textbox from virmen doesn't work. think position is wrong
% Define a textbox and set its position, size and color. 
% vr.text(1).position = [-1.2 1]; % upper-left corner of the screen
% vr.text(1).size = 0.03; % letter size as fraction of the screen
% vr.text(1).color = [1 1 0]; % yellow
% vr.text(1).window = 1; % display text in the first window


%text box works
% Define a textbox and set its position, size and color
% vr.text(1).position = [-.14 0]; % upper-left corner of the screen
% vr.text(1).size = 0.03; % letter size as fraction of the screen
% vr.text(1).color = [1 1 0]; % yellow
% vr.text(1).window = 1; % display text in the first window

%just displays "0"
% vr.text(2).string = '0';
% vr.text(2).position = [-.14 0];
% vr.text(2).size = .03;
% vr.text(2).color = [1 1 0];

%vr.friction = 0.1; % define friction that will reduce velocity by 70% during collisions

vr.startTime = now;

% --8- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)
%original single dir track code
% if (vr.position(2)>vr.endZone)
%     %vr.scaling(2) = vr.scaling(1)+(vr.scaling(2)-vr.scaling(1))*vr.scalingDecay;
%     vr.isReward = 1;
%     reward(vr,vr.timeSolenoid);
%     %%Teleport added
%     vr.position(2) = vr.beginZone;
%     vr.dp(:) = 0;
if vr.collision % test if the animal is currently in collision
    % reduce the x and y components of displacement
    %vr.dp(1:2) = vr.dp(1:2) * vr.friction;
    %vr.dp(1) = vr.dp(1) * 0.95;
    vr.dp(1) = 0;
%     vr.dp(4) = vr.dp(4)*0.2;
%     vr.position(2) = 100;
%     vr.dp(:) = 0;
end  




% % On every iteration, update the string to display the time elapsed
% vr.text(1).string = ['TIME ' datestr(now-vr.startTime,'MM.SS')];
% vr.text(1).string = num2str(vr.currentGoal);
% vr.text(2).string = ['TIME ' datestr(now-vr.startTime,'MM.SS')];
vr.text(1).string=num2str(vr.dp(4));

%from "linearTrackTwoDirections"
%symbolYPosition = 2*(vr.position(2)-vr.trackMinY)/(vr.trackMaxY-vr.trackMinY) - 1;
%vr.plot(1).y = [-1 -1 1 1 -1]*vr.symbolSize + symbolYPosition;

if (vr.currentGoal==1 && vr.position(2)>vr.topTarget) || (vr.currentGoal == 0 && vr.position(2)<vr.bottomTarget)
    vr.currentGoal = 1-vr.currentGoal;
        %vr.scaling(2) = vr.scaling(1)+(vr.scaling(2)-vr.scaling(1))*vr.scalingDecay;
     
    vr.isReward = 1;
%     reward_double(vr,vr.timeSolenoid);
% reward(vr,vr.timeSolenoid);
    vr.numRewards = vr.numRewards + 1;

    
%below for testing water reward amount
%keeps dispensing water rewards in top zone, 1Hz.
%turn off reward code above
% if (vr.position(2)>vr.topTarget)&&vr.testRew < 101;
% vr.isReward = 1;
%     reward(vr,vr.timeSolenoid);
%     vr.numRewards = vr.numRewards + 1;
%     vr.testRew = vr.testRew + 1;
%     pause(1);

         %%Teleport added
%      vr.position(2) = vr.beginZone;
%      vr.dp(:) = 0;
    
    %original again
    else
        vr.isReward = 0;
end

if vr.position(2) > vr.topTarget % test if the animal is at the end of the track (y > 200)
    %vr.dp(4) = vr.dp(4)*3.5; % set the animal’s y position to 0
    vr.dp(4) = vr.dp(4)*vr.end_gain; % set the animal’s y position to 0
    %vr.dp(:) = 0; % prevent any additional movement during teleportation
    if vr.currentGain == 2  %if in world 2, modify rot gain in track
        vr.dp(4) = vr.dp(4)*vr.RotGainFactorEnd;
    end
end

if vr.position(2) < vr.bottomTarget % test if the animal is at the end of the track (y > 200)
       %vr.dp(4) = vr.dp(4)*3.5; % set the animal’s y position to 0
    vr.dp(4) = vr.dp(4)*vr.end_gain; % set the animal’s y position to 0
    %vr.dp(:) = 0; % prevent any additional movement during teleportation
    if vr.currentGain == 2  %if in world 2, modify rot gain in track
        vr.dp(4) = vr.dp(4)*vr.RotGainFactorEnd;
    end
end
if vr.position(2) > vr.bottomTarget  && vr.position(2)<vr.topTarget
    vr.dp(4) = vr.dp(4)*vr.track_gain;
    if vr.currentGain == 2  %if in world 2, modify rot gain in track
        vr.dp(4) = vr.dp(4)*vr.RotGainFactor;
    end
    
end
    
    %teleport world
%     if double(vr.keyPressed) == 43  %ascii code for "+" %works with 2011
    if double(vr.keyPressed) == 334  %ascii code for "+" %new version with 2015 on
            vr.isReward = 1;
            
%             reward(vr,vr.timeSolenoid);
            vr.numRewards = vr.numRewards + 1;
        if vr.currentWorld == 1 
            vr.currentWorld = 2; % set the current world
        else
            vr.currentWorld = 1;
        end  
        %vr.position(2) = 0; % set the animal’s y position to 0
        if vr.currentGain == 1 
            vr.currentGain = 2; % set the current world
        else
            vr.currentGain = 1;
        end  
    end
 timestamp = now;%170512
%  % write timestamp and the x & y components of position and velocity to a file
% % using floating-point precision
% fwrite(vr.fid, [timestamp vr.position(1:2) vr.velocity(1:2)],'double');
 fwrite(vr.fid, [timestamp vr.keyPressed],'double');%170512. logs keys pressed. problem with out ofbounds. errors
%  fwrite(vr.fid, [timestamp vr.dp(4)],'double');%170512. logs keys pressed. problem with out ofbounds. errors
    


% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)
% fclose(vr.fid);%170512