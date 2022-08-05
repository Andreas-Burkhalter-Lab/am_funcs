function clo = comedibuf_read(devname, buffer, cli)
% COMEDIBUF_READ: Read streaming input (from a comedi command) into a buffer
% This function reads incoming data from the card and stores them
% in a MATLAB matrix.  Individual reads may be of undetermined
% size, but this function provides an interface for assuring that
% data may be analyzed in fixed-size chunks.
%
% Syntax:
%   clo = comedibuf_read(devname, buffer, cli)
% where
%   devname is a string, the name of the device
%   buffer is a MATLAB matrix, of size nchannels-by-nscans. It must
%     be an int16.  (Create with int16(zeros(nchannels,nscans)).)
%   cli is the input current length of the buffer; it says how many
%     scans have already been added to the buffer
% and
%   clo is the output current length of the buffer.  If the read
%     command has terminated and there are no data to return, clo
%     will be set to -1.
%
% As long as clo < prod(size(buffer)), the buffer is not yet full.
% If it's important to you to have fixed-size chunks to analyze,
% you will just want call comedibuf_read again until the buffer is
% full.  The following loop would carry this procedure out:
%   while (curlength < prod(size(buffer)))
%     curlength = comedibuf_read(devname, buffer, curlength);
%   end
%
% The buffer is full when clo == prod(size(buffer)). At this point,
%   you need to do the following:
%   (1) Reset the current length to 0, so the buffer can be filled
%       from the beginning the next time you call this function.
%   (2) Save and/or analyze the data in this buffer.
%
% If the command finishes (because you specified a fixed number
% of scans) before filling the buffer, comedibuf_read
% will resize the buffer so that prod(size(buffer)) == clo.  If
% you're going to issue a new command, you may need to reallocate
% the buffer as well, or you may be surprised by the size of the
% chunks that you receive. 
%
% See also: COMEDI_COMMAND, COMEDIBUF_WRITE, COMEDI_OPEN.
