function flsout = temporal_order(flsin)
% WHISTLES_TEMPORAL_ORDER: arrange mouse song recordings by recording time
% Syntax:
%   flsout = temporal_order(flsin)
% where
%   flsin and flsout are both cell arrays of filenames. On output, the
%     filenames are ordered so that the earliest-recorded file is first,
%     and the latest-recorded file is last.

% Copyright 2007 by Timothy E. Holy
  
  for i = 1:length(flsin)
    h = ReadAIHeaderWU1(flsin{i});
    [pathname,basename{i}] = fileparts(flsin{i});
    rectime(i) = datenum([h.date ' ' h.time],'dd-mmm-yyyy HH:MM PM');
    %fprintf('%s %s\n',h.date,h.time);
  end
  [srt,index] = sort(rectime);
  flsout = basename(index);
  