function code = Track_180cm_eh_test1_2Xrew_world_switchRotGain_CPP_test
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
vr = initializeDAQ_CPP(vr);
vr.timeSolenoid = 140; %in milliseconds
vr.topTarget = 120;%from "linearTrackTwoDirections"
vr.bottomTarget = 60;%new term
% %vr.endZone = 170;  %original
vr.beginZone = 90;

vr.rewTimeA=1880;
vr.rewTimeB=1885;

% vr.end_gain = 0.05;
vr.track_gain = 1;
vr.RotGainFactor= 1.5;%for changing rot gain in track
vr.currentGoal = 0;
%vr.scaling = [eval(vr.exper.variables.scalingGoal) eval(vr.exper.variables.scalingStart)];
%vr.topTarget = (eval(vr.exper.variables.numCylinders)-1)*eval(vr.exper.variables.cylinderSpacing);
%vr.scalingDecay = eval(vr.exper.variables.scalingDecay);
%vr.currentGoal = 1;
vr.numRewards = 0;%actual earned rewards
vr.water = 0;%keep track of total rewards for water tracking
vr.currentGain = 1;%gain world. 1=default, 2=gainFactor
vr.testRew=0;%track rew num for testing reward amount. see below
%text box
% % Define a textbox and set its position, size and color
% vr.text(1).position = [-2.5 .8]; % upper-left corner of the screen
% vr.text(1).size = 0.03; % letter size as fraction of the screen
% vr.text(1).color = [1 1 0]; % yellow

% vr.text(2).string = '0';
% vr.text(2).position = [-.14 0];
% vr.text(2).size = .03;
% vr.text(2).color = [1 1 0];

%vr.friction = 0.1; % define friction that will reduce velocity by 70% during collisions
%below for tracking time for rew.
vr.Aelapse=0;
vr.Belapse=0;
vr.Astart=0;
vr.Bstart=0;
vr.Arewards=0;
vr.Brewards=0;
%below for tracking total time
vr.Atime=0;
vr.Btime=0;
vr.curr=0;
vr.last=0;
vr.startTime = now;

% --8- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)
vr.curr=clock;
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
%     vr.dp(1) = -vr.dp(1);
    vr.dp(1) = 0;
%     vr.dp(2) = vr.dp(2)*0.1;
%     vr.dp(4) = vr.dp(4)*0.2;
%     vr.position(2) = 100;
%     vr.dp(:) = 0;  
end  

%     if vr.position(1)<= -8 %if vr.position(1)<= -trackWidth+wallWidth
%         if vr.position(2)>= 60
%             if vr.position(2)<= 120
%                 vr.position(1)=-7.9;
%             end
%         end
%     end
%     if vr.position(1)>= 8%     if vr.position(1)>= vtrackWidth-wallWidth
%         if vr.position(2)>= 60
%             if vr.position(2)<= 120
%                 vr.position(1)=7.9;
%             end
%         end
%     end

% vr.text(1).string = num2str(vr.collision);
% vr.text(1).string = num2str(vr.Aelapse);
% vr.text(1).string = num2str(vr.position(1));

if vr.position(2)<vr.bottomTarget
    vr.Atime=vr.Atime+etime(vr.curr, vr.last);
    if vr.Astart(1)>0
        vr.Aelapse=etime(vr.curr, vr.Astart);
    else
        vr.Astart=vr.curr;        
    end
end

if vr.Aelapse>vr.rewTimeA
%     disp('reward A');
    reward(vr,vr.timeSolenoid);
    vr.numRewards = vr.numRewards + 1;
    vr.Arewards=vr.Arewards+1;
%     vr.Atime=vr.Atime+vr.Aelapse;
    vr.Astart=vr.curr; 
    vr.Aelapse=0;
end

if vr.position(2)>vr.topTarget
%     vr.Bstart=0;
    vr.Btime=vr.Btime+etime(vr.curr, vr.last);
    if vr.Bstart(1)>0
        vr.Belapse=etime(vr.curr, vr.Bstart);
    else
        vr.Bstart=clock;        
    end
end

if vr.Belapse>vr.rewTimeB
%     disp('reward B');
    reward(vr,vr.timeSolenoid);
    vr.numRewards = vr.numRewards + 1;
    vr.Brewards=vr.Brewards+1;
%     vr.Btime=vr.Btime+vr.Belapse;
    vr.Bstart=vr.curr; 
    vr.Belapse=0;
end



if (vr.position(2) >vr.bottomTarget) && (vr.position(2) <vr.topTarget)
    if vr.Astart(1)>0
        vr.Aelapse=etime(vr.curr, vr.Astart);
%         vr.Atime=vr.Atime+vr.Aelapse;
        vr.Atime=vr.Atime+etime(vr.curr, vr.last);
        vr.Astart=0;
        vr.Aelapse=0;
    end
    
    if vr.Bstart(1)>0
        vr.Belapse=etime(vr.curr, vr.Bstart);
%         vr.Btime=vr.Btime+vr.Belapse;
        vr.Btime=vr.Btime+etime(vr.curr, vr.last);
        vr.Bstart=0;
        vr.Belapse=0;
    end
    
end

vr.dp(4) = vr.dp(4)*vr.track_gain;
vr.last=clock;






%text box not working
% % On every iteration, update the string to display the time elapsed
% vr.text(1).string = ['TIME ' datestr(now-vr.startTime,'MM.SS')];
% vr.text(1).string = ['TIME ' datestr(vr.startTime)];%craash
% vr.text(1).string = [num2str(vr.startTime),4];
% vr.text(1).string = [num2str(now-vr.startTime),4];%crash




%from "linearTrackTwoDirections"
%symbolYPosition = 2*(vr.position(2)-vr.trackMinY)/(vr.trackMaxY-vr.trackMinY) - 1;
%vr.plot(1).y = [-1 -1 1 1 -1]*vr.symbolSize + symbolYPosition;

% if (vr.currentGoal==1 && vr.position(2)>vr.topTarget) || (vr.currentGoal == 0 && vr.position(2)<vr.bottomTarget)
%     vr.currentGoal = 1-vr.currentGoal;
%         %vr.scaling(2) = vr.scaling(1)+(vr.scaling(2)-vr.scaling(1))*vr.scalingDecay;
%      
%     vr.isReward = 1;
%     reward_double(vr,vr.timeSolenoid);
%     vr.numRewards = vr.numRewards + 1;
% 
%     
% %below for testing water reward amount
% %keeps dispensing water rewards in top zone, 1Hz.
% %turn off reward code above
% % if (vr.position(2)>vr.topTarget)&&vr.testRew < 101;
% % vr.isReward = 1;
% %     reward(vr,vr.timeSolenoid);
% %     vr.numRewards = vr.numRewards + 1;
% %     vr.testRew = vr.testRew + 1;
% %     pause(1);
% 
%          %%Teleport added
% %      vr.position(2) = vr.beginZone;
% %      vr.dp(:) = 0;
%     
% %original again
%     else
%         vr.isReward = 0;
% end
% if vr.position(2) > vr.topTarget % test if the animal is at the end of the track (y > 200)
%     %vr.dp(4) = vr.dp(4)*3.5; % set the animal’s y position to 0
%     vr.dp(4) = vr.dp(4)*vr.end_gain; % set the animal’s y position to 0
%     %vr.dp(:) = 0; % prevent any additional movement during teleportation
% end
% 
% if vr.position(2) < vr.bottomTarget % test if the animal is at the end of the track (y > 200)
%     %vr.dp(4) = vr.dp(4)*3.5; % set the animal’s y position to 0
%     vr.dp(4) = vr.dp(4)*vr.end_gain; % set the animal’s y position to 0
%     %vr.dp(:) = 0; % prevent any additional movement during teleportation
% end
% if vr.position(2) > vr.bottomTarget  && vr.position(2)<vr.topTarget
%     vr.dp(4) = vr.dp(4)*vr.track_gain;
%     if vr.currentGain == 2  %if in world 2, modify rot gain in track
%         vr.dp(4) = vr.dp(4)*vr.RotGainFactor;
%     end
% end
%     
%     %teleport world
%     if double(vr.keyPressed) == 43  %ascii code for "+"
%             vr.isReward = 1;
%             reward(vr,vr.timeSolenoid);
%             vr.numRewards = vr.numRewards + 1;
%         if vr.currentWorld == 1; 
%             vr.currentWorld = 2; % set the current world
%         else
%             vr.currentWorld = 1;
%         end  
%         %vr.position(2) = 0; % set the animal’s y position to 0
%         if vr.currentGain == 1; 
%             vr.currentGain = 2; % set the current world
%         else
%             vr.currentGain = 1;
%         end  
%         
%         
%         
%         
%         
%     end
    
    
    
    
    
    


% if vr.isReward
%     vr.numRewards = vr.numRewards + 1;  
%     
% end





% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)
