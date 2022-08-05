function nwritten = comedibuf_write(devname, data)
% COMEDIBUF_WRITE: Write streaming data into a buffer for analog output
% This function sends data to the card for analog output.  The
% analog output command must be established using comedi_command.
%
% Syntax:
%   nwritten = comedibuf_write(devname, data)
% where
%   devname is a string, the name of the device
%   data is a MATLAB matrix, of size nchannels-by-nsamples. It must
%     be an int16.  (Create with int16(waveform).)  For a 12-bit
%     card, the values are 0 to 4095.
% and
%   nwritten is the number of samples actually written.
%
% It is not guaranteed that all the data samples will be written to
% comedi's buffer.  If a single call does not write all the
% data to comedi's buffer, you will need to keep track of the
% position within the data in a manner similar to the following:
%
%    n = 0;
%    while (n < prod(size(data)))
%      m = comedibuf_write(devname,data(1+n:end));
%      n = n+m;
%    end
%
% As long as you've set things up so the AO command doesn't terminate in the
% middle of this process, this will send the entire contents of 'data'.
%
% Note that this command writes samples to an internal buffer, and does not
% wait for them to appear on the actual output.  This function will likely
% return before all these samples have actually been output; if you
% immediately cancel the command or close the device, it may prematurely
% terminate the output of these samples.
%
% See also: COMEDI_COMMAND, COMEDIBUF_READ, COMEDI_OPEN.
