dev = '/dev/comedi0';
if ~comedi_isopen(dev)
  comedi_open(dev);   % Open the device
end
subdev = 2;
for i = 1:8
  comedi_dio_config(dev,subdev,i-1,'COMEDI_OUTPUT');
end
writemask = 255;
%comedi_dio_bitfield(dev,subdev,writemask,8);