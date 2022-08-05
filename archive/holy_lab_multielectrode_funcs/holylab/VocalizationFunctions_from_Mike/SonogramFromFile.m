function sng = SonogramFromFile(filename,nfreq,navg,plotint,imH)
% SonogramFromFile: compute the sonogram from a .bin file
% Two calling modes:
%    sng = SonogramFromFile(filename,nfreq,navg)
% produces the sonogram only. nfreq gives the number of frequencies to use,
% and navg is the number of consecutive blocks to average together.
%    sng = SonogramFromFile(filename,nfreq,navg,plotint,imH)
% plots the sonogram as it is computed. The plotting interval is
% determined by plotint, and imH contains the handle to the image
% object, which must already be allocated at the appropriate size.
%
% See also SFF.
