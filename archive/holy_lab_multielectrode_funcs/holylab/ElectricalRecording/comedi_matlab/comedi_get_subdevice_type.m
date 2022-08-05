function type = comedi_get_subdevice_type(devname,subdevice)
% COMEDI_GET_SUBDEVICE_TYPE: find the appropriate comedi subdevice
% Syntax:
%    type = comedi_get_subdevice_type(devname,subdevice)
% where
%   devname is a string, the name of the device
%   subdevice is an integer, the number of the subdevice
% and
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
%
% See also: COMEDI_OPEN, COMEDI_FIND_SUBDEVICE_BY_TYPE.
