function write_fake_merec(fileoutname,w,headerfilename)
% WRITE_FAKE_MEREC: create a .merec file from voltage traces
% Syntax:
%   write_fake_merec(fileoutname,w,headerfilename)
% where
%   fileoutname is the name of the file (usually with .merec extension)
%     that you want to save the waveform data to
%   w is a n_channels-by-n_scans matrix containing the waveform data
%   headerfilename is the name of a text file that contains a working
%     example of the merec header format (snipped out from the beginning
%     of a real .merec file)
%
% On output, the file will contain data padded up to the number of
% channels listed in the .merec header's channel list. The "digitization"
% parameters will be set to optimally capture the range of voltages in
% w. The comments field will be populated to indicate that the file is
% written by write_fake_merec.
%
% See also: SYNTHESIZE_WAVEFORM.
  
% Copyright 2007 by Timothy E. Holy
  
  [fid,msg] = fopen(headerfilename,'r');
  if (fid < 0)
    error(['Error opening ' headerfilename ': ' msg])
  end
  headertxt = fread(fid,Inf,'char');
  fclose(fid);
  headertxt = char(headertxt');

  % Check # of channels, and pad w if necessary
  channel_list = str2num(key2value(headertxt,'channel list'));
  n_channels = length(channel_list);
  if (size(w,1) < n_channels)
    warning('ElectricalRecording:FewerChannels',...
      'Added zeros to bring # of channels from %d to %d\n',...
      size(w,1),n_channels);
    w(n_channels,1) = 0;
  end

  % Create an 12-bit integer (or whatever) version of the waveform,
  % spanning the range of recorded values
  w_min = min(w(:));
  w_max = max(w(:));
  min_sample = str2num(key2value(headertxt,'min sample'));
  max_sample = str2num(key2value(headertxt,'max sample'));
  m = (max_sample-min_sample)/(w_max-w_min);
  w_scaled = m*(w-w_min) + min_sample;
  
  n_scans = size(w,2);
  % Adjust several components of the header text
  [h1,h2] = split_header(headertxt,'min input');
  headertxt = [h1 num2str(w_min) h2];
  [h1,h2] = split_header(headertxt,'max input');
  headertxt = [h1 num2str(w_max) h2];
  [h1,h2] = split_header(headertxt,'nscans');
  headertxt = [h1 num2str(n_scans) h2];
  [h1,h2] = split_header(headertxt,'stimulus sequence');
  headertxt = [h1 h2];
  [h1,h2] = split_header(headertxt,'datetime');
  headertxt = [h1 datestr(now,'MM/DD/YY ') datestr(now,'HH:MM:SS PM') h2];
  [h1,h2] = split_header(headertxt,'comment');
  headertxt = [h1 'Written by write_fake_merec' h2];
  headersz = numel(headertxt)+10;
  [h1,h2] = split_header(headertxt,'header size');
  headertxt = [h1 num2str(headersz) h2];
  if (numel(headertxt) >= headersz)
    error('Need more padding on header size');
  end
  headertxt(headersz) = 0;  % insert padding to make correct size
  
  [fid,msg] = fopen(fileoutname,'w');
  if (fid < 0)
    error(['Error opening ' fileoutname ' for writing: ' msg])
  end
  count = fwrite(fid,headertxt,'char');
  if (count ~= headersz)
    error('Error writing header');
  end
  if (min_sample == 0)
    count = fwrite(fid,w_scaled,'uint16');
  else
    count = fwrite(fid,w_scaled,'int16');
  end
  fclose(fid);
  
function [h1,h2] = split_header(txt,key)
  tExpression= strcat('^', key,'=([^\n]*)$');
  [s, f, t]=regexp(txt, tExpression, 'once', 'lineanchors');
  h1 = txt(1:t(1)-1);
  h2 = txt(f+1:end);
  