ochan = 0;
ichan = 0;
wavefreq = 1000;
oscanrate = 50000;
iscanrate = 50000;
intime = 1; % Input time in seconds
amplitude = 1; % Output amplitude, in volts
% Create the sinusoidal waveform
% Let the waveform last 10% longer than the intime, so there is
% no starting and stopping
extrafac=1.1;
npts = ceil(extrafac*intime*oscanrate);
x = linspace(0,extrafac*intime*wavefreq*2*pi,npts);
y = int16(amplitude/10*2048*sin(x) + 2048);  % Card (on range 0 ) goes from -10V to +10V
devname = '/dev/comedi0';
osubdev = 1;
isubdev = 0;
orange = 0;  % If you change this, you'll have to fix some other code here
irange = 0;
toclose = 0;
if ~comedi_isopen(devname)
    comedi_open(devname);
    toclose = 1;
end
%insn = [];
%insn(1).insn = 'INSN_INTTRIG';
%insn(1).insn = 'INSN_WAIT';
%insn(1).data = 9998;
%insn(1).n = 1;
%insn(1).subdev = 0;
%insn(2).insn = 'INSN_INTTRIG';
%insn(2).subdev = 0;
%insn(3).insn = 'INSN_INTTRIG';
%insn(3).subdev = 1;
cmdao = cmd_ao_trig(osubdev,ochan,oscanrate,orange);
%cmdao.stop_arg = length(y);
%cmdao.stop_src = 'TRIG_NONE';
%cmdao.stop_arg = 0;
ibuf = int16(zeros(1,intime*iscanrate));
cmdai = cmd_ai_perpetual(isubdev,ichan,iscanrate,irange);
cmdai.stop_src = 'TRIG_COUNT';
cmdai.stop_arg = length(ibuf);
cmdai.start_src = 'TRIG_INT';
msg = 'failure';
while ~strcmp(msg,'success')
  [msg,cmdai] = comedi_command_test(devname,cmdai);
end
%'issue ai command'
if (comedi_command(devname,cmdai) < 0)
    error('input command failed');
end
if (comedi_command(devname,cmdao) < 0)
    error('output command failed');
end
% Load output buffer with some ready-to-go samples
nwritten = comedibuf_write(devname,y);
% Trigger the commands
%ilout = comedi_do_insnlist(devname,insn)
if (comedi_internal_trigger(devname,osubdev,0) < 0)
    error('Error triggering output');
end
if (comedi_internal_trigger(devname,isubdev,0) < 0)
    error('Error triggering input');
end
% Read & write simultaneously
cl = 0;
while (cl < prod(size(ibuf)))
  cl = comedibuf_read(devname,ibuf,cl);
  while (nwritten < prod(size(y)))
      m = comedibuf_write(devname,y(1+nwritten:end));
      nwritten = nwritten+m;
  end
end
% Stop the analog output command
comedi_cancel_and_flush(devname,osubdev);
% Set analog output to zero volts
comedi_data_write(devname,osubdev,ochan,0,'AREF_GROUND',2048);
% Plot the data
[range,unit] = comedi_get_range(devname,isubdev,ichan,irange);
ym = diff(range)*double(ibuf)/4095 + range(1);
plot(linspace(0,intime,length(ym)),ym);
% Close device if needed
if toclose
    comedi_close(devname)
end