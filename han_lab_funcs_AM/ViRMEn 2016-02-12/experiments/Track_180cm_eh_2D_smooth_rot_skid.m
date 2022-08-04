function code = Track_180cm_eh_2D_smooth_rot_skid

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
% vr = initializeDAQ(vr);
vr = initializeDAQ_smooth_EH(vr);%for smoothing ball input. used with moveBall_2D_smooth
vr.timeSolenoid = 140; %in milliseconds
vr.topTarget = 170;%from "linearTrackTwoDirections"
vr.bottomTarget = 10;%new term
%vr.endZone = 170;  %original
vr.beginZone = 10;

vr.currentGoal = 0;
%vr.scaling = [eval(vr.exper.variables.scalingGoal) eval(vr.exper.variables.scalingStart)];
%vr.topTarget = (eval(vr.exper.variables.numCylinders)-1)*eval(vr.exper.variables.cylinderSpacing);
%vr.scalingDecay = eval(vr.exper.variables.scalingDecay);
%vr.currentGoal = 1;
vr.numRewards = 0;

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
    vr.dp(4)=vr.dp(4)*.5;
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

if (vr.currentGoal==1 && vr.position(2)>vr.topTarget) || (vr.currentGoal == 0 && vr.position(2)<vr.bottomTarget)
    vr.currentGoal = 1-vr.currentGoal;
        %vr.scaling(2) = vr.scaling(1)+(vr.scaling(2)-vr.scaling(1))*vr.scalingDecay;
     
    vr.isReward = 1;
    reward(vr,vr.timeSolenoid);
    vr.numRewards = vr.numRewards + 1;
    
    
%below for testing water reward amount
%keeps dispensing water rewards in top zone, 1Hz
% if (vr.position(2)>vr.topTarget);
% vr.isReward = 1;
%     reward(vr,vr.timeSolenoid);
%     vr.numRewards = vr.numRewards + 1;
%     pause(1);

         %%Teleport added
%      vr.position(2) = vr.beginZone;
%      vr.dp(:) = 0;
    
%original again
else
    vr.isReward = 0;
if vr.position(2) > 175 % test if the animal is at the end of the track (y > 200)
    vr.dp(4) = vr.dp(4)*3; % set the animal�s y position to 0
%     vr.dp(4) = vr.dp(4)*2; % set the animal�s y position to 0
    %vr.dp(:) = 0; % prevent any additional movement during teleportation
end

if vr.position(2) < 5 % test if the animal is at the end of the track (y > 200)
    vr.dp(4) = vr.dp(4)*3; % set the animal�s y position to 0
%     vr.dp(4) = vr.dp(4)*2; % set the animal�s y position to 0
    %vr.dp(:) = 0; % prevent any additional movement during teleportation
end
    
    %teleport world
    if double(vr.keyPressed) == 43  %ascii code for "+"
            vr.isReward = 1;
            reward(vr,vr.timeSolenoid);
            vr.numRewards = vr.numRewards + 1;
        if vr.currentWorld == 1; 
            vr.currentWorld = 2; % set the current world
        else
            vr.currentWorld = 1;
        end  
        %vr.position(2) = 0; % set the animal�s y position to 0
        
    end
    
    
    
    
    
    
end

% if vr.isReward
%     vr.numRewards = vr.numRewards + 1;  
%     
% end





% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)
