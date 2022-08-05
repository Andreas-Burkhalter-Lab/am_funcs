% An example script showing streaming analog output
devname = '/dev/comedi0';
sampfreq = 1000;   % Output 1000 pts/sec
wavefreq = 10;     % A sine wave of 10Hz
subdevice = 1;  % This has to be an analog output device. See aitest for
                % examples of how to find such a device within MATLAB.
channels = 0;   % Can use more than one channel, as long as you
                % construct appropriate waveforms
outtime = 10;   % output samples for at least 10s
% Open device
toclose = 0;
if ~comedi_isopen(devname)
    comedi_open(devname);
    toclose = 1;
end
% Prepare output command to run perpetually;
% will require a trigger to launch it
cmd = cmd_ao_trig(subdevice,channels,sampfreq,0);
% Queue the command
[msg,cmdo] = comedi_command_test(devname,cmd);
cmd=cmdo;
cmdval = comedi_command(devname,cmd);
if (cmdval < 0)
  error('Command failed!');
end
% Prepare a waveform. We might make several wavelengths worth,
% so data can be sent in bigger chunks to the buffer.
wavelen = round(sampfreq/wavefreq);
minsamples = 0.2*sampfreq + 2048;  % Gives a minimum size, plus requires update no more frequently than 0.2s
nper = ceil(minsamples/wavelen);
nper, wavelen
x = linspace(0,nper*2*pi,nper*wavelen);  % nper wavelengths
y = int16(2048+sin(x)*204);    % Assumes at least 12 bit AO
% Send the data
m = comedibuf_write(devname,y);
fprintf('Pre-send %d samples\n',m);
fprintf('Writing %ds worth of samples...',outtime);
% Trigger the command
if (comedi_internal_trigger(devname,subdevice,0) < 0)
  error('Trigger failed');
end
% Keep sending data
tic
starttime = toc;
n = m;
while (toc < starttime+outtime)
  if (n < prod(size(y)))
    m = comedibuf_write(devname,y(1+n:end));
    n = n+m;
  else
    n = 0;
  end
end
comedi_cancel_and_flush(devname,subdevice);
% Restore output value to 0V
comedi_data_write(devname,subdevice,channels,0,'AREF_GROUND',0*channels + 2048);
if (toclose)
    comedi_close(devname)
end
