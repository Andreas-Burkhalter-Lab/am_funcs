function subdev = comedi_find_subdevice_by_type(devname,type,start_subdevice)
% COMEDI_FIND_SUBDEVICE_BY_TYPE: find the appropriate comedi subdevice
% Syntax:
%    subdev = comedi_find_subdevice_by_type(devname,type,start_subdevice)
% where
%   devname is a string, the name of the device
%   type is a string which may take the following values:
%     COMEDI_SUBD_AI (analog input)
%     COMEDI_SUBD_AO (analog output)
%     COMEDI_SUBD_DI (digital input)
%     COMEDI_SUBD_DO (digital output)
%     COMEDI_SUBD_DIO (digital input/output)
%     COMEDI_SUBD_COUNTER (counter)
%     COMEDI_SUBD_TIMER (timer)
%     COMEDI_SUBD_MEMORY (memory, EEPROM, DPRAM)
%     COMEDI_SUBD_CALIB (calibration DACs)
%     COMEDI_SUBD_PROC (processor, DSP
%   start_subdevice (optional) determines the starting subdevice
%     number, the search proceding incrementally.  If absent, this
%     defaults to 0.
% and
%   subdev is the subdevice number.
%
% This function gives an error if a subdevice of the requested type
% cannot be found.
%
% See also: COMEDI_OPEN, COMEDI_GET_SUBDEVICE_TYPE.
