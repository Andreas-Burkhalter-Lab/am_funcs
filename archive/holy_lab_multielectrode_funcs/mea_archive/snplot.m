%%%% for looking at spikes and motion noise in raw recording
close all
figure
% merecfile = '3sizetuning.merec';
% merecfile = '1rf.merec';
merecfile = 'sftfss.merec';
m = merecmm(merecfile);
trange = [1000 1010]; % secs
chan = 34;
scanrange = [1+m.scanrate*trange(1) : m.scanrate*trange(2)];
trace = m(chan,scanrange);
timeaxis = scanrange / m.scanrate; 
plot(timeaxis,trace);
