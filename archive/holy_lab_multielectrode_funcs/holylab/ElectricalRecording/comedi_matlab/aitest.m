devname = '/dev/comedi0';
% Use the analog input subdevice
subdevtype = 'COMEDI_SUBD_AI';
readchannels = 0:3;   % Read from channels 0 to 3, inclusive.
plotchannels = 0:1;   % Plot channels 0 (blue) and 1 (green).
separateplots = 0;
freq = 10000;      % Sampling frequency (in scans/sec)
nbufs_tot = 1000;
screentime = 1;    % Number of times per second to draw on the screen
if ~comedi_isopen(devname)
  comedi_open(devname);   % Open the device
end
% Find the subdevice number from the type
subdev = comedi_find_subdevice_by_type(devname,subdevtype);
% Check to make sure that we don't ask for too many channels
maxchannels = comedi_get_n_channels(devname,subdev);
if (max(readchannels)+1 > maxchannels)
  readchannels = 0:(maxchannels-1);
  warning(['Only ' num2str(maxchannels) ' supported on subdevice']);
end
% Compute the index for channel plotting. The approach here
% insures the color order will respect the order in plotchannels,
% even if not monotonic
[pc,ipc,irc] = intersect(plotchannels,readchannels);
[sipc,permindx] = sort(ipc);
plotindex = irc(permindx);
pc = pc(permindx);
% Increase the buffer size to the maximum
bufsz = comedi_get_max_buffer_size(devname,subdev);
comedi_set_buffer_size(devname,subdev,bufsz);
bufsz
% Try to record from -5 to +5 volts, or as near as card will allow
rangeindex = comedi_find_range(devname,subdev,0,'UNIT_volt',[-5 5]);
range = comedi_get_range(devname,subdev,0,rangeindex);
% Set up comedi command
cmd = cmd_ai_perpetual(subdev,readchannels,freq,rangeindex);
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
% Prepare the read buffer
readbuf = int16(zeros(length(readchannels),freq*screentime));
% Prepare the plot
hfig = figure;
if separateplots
  for i = 1:length(plotindex)
    hax(i) = subplot(length(plotindex),1,i);
  end
else
  hax = gca;
end
set(hax,'YLim',range,'XLim',[1 size(readbuf,2)],'NextPlot','replacechildren');
ylabel('Voltage (V)');
xlabel('Sample number');
drawnow     % open the plot window to get this time-consuming step out
            % of the way

% Execute the command
% After this, the board is actually running, and you need to empty the
% read buffer sufficiently quickly
if (comedi_command(devname,cmd) < 0)
  error(['Command failed. Perhaps there is a command already running?', ...
         '  Try comedi_cancel_and_flush(devname,subdev)']);
end
% Now start reading
nbufs = 0;
curlength = 0;
while (nbufs < nbufs_tot & ishandle(hax))
  curlength = comedibuf_read(devname,readbuf,curlength);
  if (curlength < 0)
    error('Read command finished by itself!');
  end
  if (curlength == prod(size(readbuf)))
    % Buffer is full, do something with it
    % Convert to volts
    plotbuf = double(readbuf) * diff(range)/4095 + range(1);  % This
                                                              % assumes a
                                                              % 12-bit card!
    % Plot it
    if (length(hax) == 1)
      plot(plotbuf(plotindex,:)','Parent',hax);
      title(['Buffer ' num2str(nbufs+1) ' out of ' ...
             num2str(nbufs_tot)],'Parent',hax);
      legend(num2str(pc'))
    else
      for i = 1:length(hax)
        plot(plotbuf(plotindex(i),:),'Parent',hax(i))
        title(['Channel ' num2str(pc(i)) ', buffer ' num2str(nbufs+1) ...
               ' out of ' num2str(nbufs_tot)],'Parent',hax(i));
      end
    end
    nbufs = nbufs + 1;
    drawnow;
    curlength = 0;
  end
end
% Terminate the command
comedi_cancel_and_flush(devname,subdev);
comedi_close(devname);
