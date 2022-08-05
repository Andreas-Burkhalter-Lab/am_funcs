% Data acquisition functions (comedi interface)
% Comedi is a general interface for many different data acquisition
% cards under Linux.  See http://stm.lbl.gov/comedi.
% These functions allow you to access some of the comedi library
% from MATLAB.
%
% Device management
%   comedi_open          - Open device 
%   comedi_isopen        - Query whether device is open
%   comedi_close         - Close device
%   comedi_find_subdevice_by_type  - Find subdevice numbers
%   comedi_get_subdevice_type      - Find type given a subdevice number
%   comedi_get_max_buffer_size     - Useful for streaming data
%   comedi_set_buffer_size         - Adjust the buffer size (streaming)
%   comedi_get_n_ranges  - Determine number of A/D voltage ranges
%   comedi_get_range     - Determine the nth A/D voltage range
%   comedi_find_range    - Find corresponding range index
%   comedi_get_n_channels - Determine the number of channels in subdevice
%   comedi_get_maxddata    - get the max value in digital unit
%
% Streaming data (for input & output, depending on hardware support)
%   comedi_command       - Initiate a board-autonomous command (streaming)
%   comedi_command_test  - Check and adjust parameters of a command
%   comedi_cancel_and_flush        - Terminate a command
%
% Instructions (useful for triggers and the like)
%   comedi_do_insnlist   - Execute an instruction list
%   comedi_internal_trigger  - Start a waiting command running
%
% Data acquisition
%   comedi_data_read     - Perform a single A/D conversion
%   comedibuf_read       - Acquire data from a command
%
% Analog output
%   comedi_data_write    - Perform a single D/A conversion
%   comedibuf_write      - Send streaming data to output command process
%
% Digital sample control
%   comedi_dio_config    - Configure digital I/O channel direction
%   comedi_dio_bitfield  - Set and get multiple digital channels simultaneously
%
% Example MATLAB scripts
%   aitest               - Analog input example using commands
%   aotest               - Analog output example using commands
%   aiaotest             - Simultaneous analog input/output example
%   diotest              - Digital I/O example
%   cmd_ai_perpetual     - Command preparation for AI
%   cmd_ao_once          - Command preparation for AO
%
% Other
%   comedi_missing       - Only part of the functionality of comedilib has been wrapped for
  %MATLAB.  This file documents some of these shortcomings and how
  %get around some of them.

