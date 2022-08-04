function reward(vr,rewardTime)
queueOutputData(vr.waterSession,[ones(1,rewardTime)*vr.highVoltage vr.lowVoltage]');
%startForeground(vr.waterSession);%runs data acq in foreground. pause during rewards
startBackground(vr.waterSession);%no pause