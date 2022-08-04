function code = Track_180cm_eh_test1_2Xrew_world_switch_noVR_randRew
%switches to black world
%normal VR rewrd off, random reward
%gaussian dist or intervals, 
% can have 1X, 2x, or 3x int. between rew
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
vr.rand1 = (5 + (15 - 5)*rand(500,1)); %5 to 15, 500 numbers
vr.rand1(randi(500,50,1))= 2*(5 + (15 - 5)*rand(50,1));
vr.rand1(randi(500,10,1))= 3*(5 + (15 - 5)*rand(10,1));

% vr.rand2 = rand(500,1);
vr = initializeDAQ(vr);
vr.timeSolenoid = 140; %in milliseconds
vr.topTarget = 170;%from "linearTrackTwoDirections"
vr.bottomTarget = 10;%new term
% %vr.endZone = 170;  %original
vr.beginZone = 8.5;
vr.end_gain = 3;
vr.track_gain = .4;
vr.RotGainFactor= 2;%for changing rot gain in track
vr.currentGoal = 0;
%vr.scaling = [eval(vr.exper.variables.scalingGoal) eval(vr.exper.variables.scalingStart)];
%vr.topTarget = (eval(vr.exper.variables.numCylinders)-1)*eval(vr.exper.variables.cylinderSpacing);
%vr.scalingDecay = eval(vr.exper.variables.scalingDecay);
%vr.currentGoal = 1;
vr.numRewards = 0;%actual earned rewards
vr.water = 0;%keep track of total rewards for water tracking
vr.currentGain = 1;%gain world. 1=default, 2=gainFactor
vr.testRew=0;%track rew num for testing reward amount. see below
%text box not working
% % Define a textbox and set its position, size and color
% vr.text(1).position = [-1.2 1]; % upper-left corner of the screen
% vr.text(1).size = 0.03; % letter size as fraction of the screen
% vr.text(1).color = [1 1 0]; % yellow

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



%text box not working
% % On every iteration, update the string to display the time elapsed
% vr.text(1).string = ['TIME ' datestr(now-vr.startTime,'MM.SS')];

% vr.text(2).string = ['TIME ' datestr(now-vr.startTime,'MM.SS')];


%from "linearTrackTwoDirections"
%symbolYPosition = 2*(vr.position(2)-vr.trackMinY)/(vr.trackMaxY-vr.trackMinY) - 1;
%vr.plot(1).y = [-1 -1 1 1 -1]*vr.symbolSize + symbolYPosition;
if vr.currentWorld == 1;  %shut off rew in vr.world =2
    if (vr.currentGoal==1 && vr.position(2)>vr.topTarget) || (vr.currentGoal == 0 && vr.position(2)<vr.bottomTarget)
        vr.currentGoal = 1-vr.currentGoal;
            %vr.scaling(2) = vr.scaling(1)+(vr.scaling(2)-vr.scaling(1))*vr.scalingDecay;
        vr.isReward = 1;
        reward_double(vr,vr.timeSolenoid);
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
end

if vr.position(2) > vr.topTarget % test if the animal is at the end of the track (y > 200)
    %vr.dp(4) = vr.dp(4)*3.5; % set the animal’s y position to 0
    vr.dp(4) = vr.dp(4)*vr.end_gain; % set the animal’s y position to 0
    %vr.dp(:) = 0; % prevent any additional movement during teleportation
end

if vr.position(2) < vr.bottomTarget % test if the animal is at the end of the track (y > 200)
    %vr.dp(4) = vr.dp(4)*3.5; % set the animal’s y position to 0
    vr.dp(4) = vr.dp(4)*vr.end_gain; % set the animal’s y position to 0
    %vr.dp(:) = 0; % prevent any additional movement during teleportation
end
if vr.position(2) > vr.bottomTarget  && vr.position(2)<vr.topTarget
    vr.dp(4) = vr.dp(4)*vr.track_gain;
    if vr.currentGain == 2  %if in world 2, modify rot gain in track
        vr.dp(4) = vr.dp(4)*vr.RotGainFactor;
    end
end

if vr.currentWorld == 2;
    vr.rewInt = vr.rand1(vr.numRewards);
    vr.curr=clock;
    if etime(vr.curr,vr.Start) > vr.rewInt
        vr.isReward = 1;
        reward_double(vr,vr.timeSolenoid);
        vr.numRewards = vr.numRewards + 1;
        vr.Start=clock;
        else
            vr.isReward = 0;
    end
    
end

    %teleport world
    if double(vr.keyPressed) == 43  %ascii code for "+"
            vr.isReward = 1;
            reward(vr,vr.timeSolenoid);
            vr.numRewards = vr.numRewards + 1;
        if vr.currentWorld == 1; 
            vr.currentWorld = 2; % set the current world
            vr.worlds{vr.currentWorld}.surface.visible(:) = false;
            vr.Start=clock;
        else
            vr.currentWorld = 1;
            vr.worlds{vr.currentWorld}.surface.visible(:) = true;
        end  
        %vr.position(2) = 0; % set the animal’s y position to 0
        if vr.currentGain == 1; 
            vr.currentGain = 2; % set the current world
        else
            vr.currentGain = 1;
        end  
    end
    
%testing ball speed mod
% vr.dp(4) = 0;%set rotation dp =0
%below crashes because of reward
% if vr.position(2) > 176 % test if the animal is at the end of the track (y > 200)
%     vr.position(2) = 4; % set the animal’s y position to 0
%     vr.dp(:) = 0; % prevent any additional movement during teleportation
% end

% if vr.isReward
%     vr.numRewards = vr.numRewards + 1;  
%     
% end





% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)
