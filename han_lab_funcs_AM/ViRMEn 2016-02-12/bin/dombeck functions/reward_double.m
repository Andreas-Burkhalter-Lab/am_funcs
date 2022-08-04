function reward_double(vr,rewardTime)
%queueOutputData(vr.waterSession,[ones(1,rewardTime)*vr.highVoltage vr.lowVoltage]');
rew_2X=ones(1,(rewardTime*2)+501);
rew_2X(1,rewardTime+1:rewardTime+501)=0;
queueOutputData(vr.waterSession,[rew_2X*vr.highVoltage vr.lowVoltage]');
%startForeground(vr.waterSession);%runs data acq in foreground. pause during rewards
startBackground(vr.waterSession);%no pause