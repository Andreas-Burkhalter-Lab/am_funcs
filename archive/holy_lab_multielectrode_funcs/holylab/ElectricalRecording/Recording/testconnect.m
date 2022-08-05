% A script to check that all channels are connected to the A/D
% card.  Goes through all channels and records some data.  Because
% disconnected channels don't drain the charge from the previous
% channel very effectively (and input is high-impedance), channels
% are read one-at-a-time.  Disconnected channels still show signal
% (residual conductance through the MUX?), but the amplitude is
% noticably lower.
%
% It's also possible to do such tests with a device with lower
% input impedance, e.g. an oscilloscope.  Still, there are virtues
% in having a quick, semi-automated check available.
nscans = 1000;
devname = '/dev/comedi0';
comedi_open(devname);
for i = 1:64
  for j = 1:nscans
    data(i,j) = comedi_data_read(devname,0,i-1,0,'AREF_GROUND');
  end
end
figure
hax = zeros(1,64);
for i = 1:64
  hax(i) = subplot(8,8,i);
  plot(data(i,:));
  title(num2str(i-1))
end
set(hax,'YLim',[min(min(data)) max(max(data))]);
comedi_close(devname);
return
