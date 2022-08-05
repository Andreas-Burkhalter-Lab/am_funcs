function flsout = sort_ai_by_time(flsin,options)
% SORT_AI_BY_TIME: order .merec, .ssnp, etc files by recording time
% Syntax:
%   flsout = sort_ai_by_time(flsin)
%   flsout = sort_ai_by_time(flsin,options)
% where
%   flsin is a cell array of file names;
%   options can be passed to readheader (e.g., for legacy file support);
% and
%   flsout is a cell array of the same file names, but sorted in
%     increasing temporal order.
%
% See also: DIRBYTIME.
  
  for i = 1:length(flsin)
    if (nargin > 1)
      h = readheader(flsin{i},options);
    else
      h = readheader(flsin{i});
    end
    dn(i) = datenum([h.date ' ' h.time]);
  end
  [sdn,sortIndex] = sort(dn);
  flsout = flsin(sortIndex);
  