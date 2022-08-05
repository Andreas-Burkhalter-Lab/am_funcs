function bitsout = comedi_dio_bitfield(devname,subdev,write_mask,bitsin)
% COMEDI_DIO_BITFIELD: set and get values from digital channels
% Syntax:
%   bitsout = comedi_dio_bitfield(devname,subdev,write_mask,bitsin)
% where
%   devname is a string, the name of the device
%   subdev is the subdevice number (see COMEDI_FIND_SUBDEVICE_BY_TYPE)
%   write_mask is an integer; for each bit equal to 1, the
%     corresponding digital channel (least significant bit =
%     channel 0) will be changed to reflect the bit value in
%     bitsin.  Bits set to 0 will result in no change of output in
%     the corresponding digital channel. 
%   bitsin is an integers; the bit values specify the
%     state of the corresponding digital channels, subject to the
%     masking from write_mask
% and
%   bitsout is the value read from all digital channels, after
%     completion of the writing step.  Channels configured for
%     output (see COMEDI_DIO_CONFIG) yield an undefined result for
%     the read operation.
%
% See also: COMEDI_OPEN, COMEDI_DIO_CONFIG, COMEDI_FIND_SUBDEVICE_BY_TYPE.
