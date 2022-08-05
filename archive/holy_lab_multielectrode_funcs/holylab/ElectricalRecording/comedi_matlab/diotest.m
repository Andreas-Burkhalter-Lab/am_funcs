devname = '/dev/comedi0';
% Use the digital I/O subdevice
subdevtype = 'COMEDI_SUBD_DIO';
if ~comedi_isopen(devname)
  comedi_open(devname);   % Open the device
end
% Find the subdevice number and the number of channels
subdev = comedi_find_subdevice_by_type(devname,subdevtype);
nchannels = comedi_get_n_channels(devname,subdev);
% Configure DIO channels for output and set them to alternate state
writemask = 0;
bits = 0;
for i = 1:nchannels
  comedi_dio_config(devname,subdev,i-1,'COMEDI_OUTPUT');
  if (mod(i,2))
    bits = bitset(bits,i,1);
  end
  writemask = bitset(writemask,i,1);
end
% Write the values to digital channels
comedi_dio_bitfield(devname,subdev,writemask,bits);
disp('Digital channels should now have alternating state.');
comedi_close(devname);
