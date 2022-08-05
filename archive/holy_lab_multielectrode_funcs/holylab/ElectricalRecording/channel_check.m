devname = '/dev/comedi0';
subdev = 0;
freq = 10000;
nsamp = 1000;
chans = 0:63;
if ~comedi_isopen(devname)
  comedi_open(devname);   % Open the device
end
bufsz = comedi_get_max_buffer_size(devname,subdev);
comedi_set_buffer_size(devname,subdev,bufsz);
% Try to record from -5 to +5 volts, or as near as card will allow
rangeindex = comedi_find_range(devname,subdev,0,'UNIT_volt',[-5 5]);
range = comedi_get_range(devname,subdev,0,rangeindex);
% Set up comedi command
cmd = cmd_ai_perpetual(subdev,0,freq,rangeindex);
cmd.stop_src = 'TRIG_COUNT';  % Terminate after nsamp samples
cmd.stop_arg = nsamp;
% Tweak command to comply with values actually supported by the board
msg = 'fail';
nfailed = 0;
maxfail = 5;
while (~strcmp(msg,'success') & nfailed < maxfail)
    [msg,cmd] = comedi_command_test(devname,cmd);  % Tunes parameters to realizable
    nfailed = nfailed + 1;
end
if (nfailed == maxfail)
    error('Command preparation did not succeed');
end
% Recreate the command so the channel can be changed
cmd = rmfield(cmd,'packedchanlist');
cmd.rangelist = rangeindex;
cmd.areflist = 0;
% Prepare the read buffer
readbuf = int16(zeros(1,nsamp));
% Loop over channels
for i = 1:length(chans)
    % Execute the command
    % After this, the board is actually running, and you need to empty the
    % read buffer sufficiently quickly
    cmdchan = cmd;
    cmdchan.chanlist = chans(i);
    if (comedi_command(devname,cmdchan) < 0)
        error(['Command failed. Perhaps there is a command already running?', ...
         '  Try comedi_cancel_and_flush(devname,subdev)']);
    end
    % Now read the data
    curlength = 0;
    while (curlength < prod(size(readbuf)))
        curlength = comedibuf_read(devname,readbuf,curlength);
        if (curlength < 0)
            error('Read command finished by itself!');
        end
    end
    % Fix for comedi error
    %wave = double(readbuf);
    %neg = find(wave < 0);
    %wave(neg) = double(readbuf(neg)) + 4096;
    waveform(:,i) = double(readbuf') * diff(range)/4095 + range(1);
end
% Plot the data
for i = 1:length(chans)
    subplot(8,8,i);
    plot(waveform(:,i));
    set(gca,'YLim',range);
    title(num2str(chans(i)));
end
comedi_close(devname);
