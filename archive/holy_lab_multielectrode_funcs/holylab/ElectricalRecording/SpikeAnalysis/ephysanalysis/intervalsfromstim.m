function [intervals,identity,vlvsout] = intervalsfromstim(ephysin,trange,vlvs)
% [intervals,identity,vlvsout] = intervalsfromstim(ephysin,trange,vlvs)
% where
%   ephysin is a structure array of type EPHYS
%   trange is in seconds, and gives the [tbegin tend] relative to valve openings
%     Or, if trange is set to the string 'open', just the time range in which the
%     valve is open is used.
%   vlvs (optional) is the list of good valve numbers (if absent or set to 'all',
%     chooses all non-zero values).
% and
%   intervals is a cell array (of length(ephysin)) of n-by-2 matrices, giving
%     start and end scan numbers of periods surrounding valve openings (periods
%     defined by trange). The intervals are in increasing temporal order.
%   identity is a cell array of the same size, each holding a vector of the
%     valve number associated with each interval
%   vlvsout is a vector of unique valve numbers that were used
%
% See also: EPHYS, EPHYSSUBRANGE.
  
% Changelog:
%   TEH, 05-03-2004: fetch stimulus info if needed
  if ~isstruct(ephysin)
    error('ephysin must be a structure (see EPHYS)');
  end
  if ~isfield(ephysin,'stimulus')
    if isfield(ephysin,'stimulusfile')
      ephysin = ephysfetch(ephysin,'stimulus');
    else
      error(['Can''t calculate intervals because stimulus information is ' ...
             'unavailable.']);
    end
  end
  nfiles = length(ephysin);
  % Find the number of valves
  if (nargin > 2 & isnumeric(vlvs))
    nvalves = length(vlvs);
    vlvsout = vlvs;
  elseif (nargin < 3 | (ischar(vlvs) & strcmp(vlvs,'all')))
    allstim = cat(2,ephysin.stimulus);
    vlvs = unique(allstim(1,:));
    vlvs = setdiff(vlvs,0);                % Do not include 0
    nvalves = length(vlvs);
    vlvsout = vlvs;
  end
  % Check input format on trange
  justopen = 0;
  if (ischar(trange) & strcmp(trange,'open'))
    justopen = 1;
  end
  identity = cell(1,nfiles);
  for i = 1:nfiles
    ikeep = cell(1,nvalves);
    for j = 1:nvalves
      ikeep{j} = find(ephysin(i).stimulus(1,:) == vlvs(j));
    end
    ikeep = sort([ikeep{:}]); % put in temporal order
    if ~justopen
      trange_scans = trange*ephysin(i).scanrate;
      temp = repmat(ephysin(i).stimulus(2,ikeep)',1,2) + ...
        repmat(trange_scans,length(ikeep),1);
      identity{i} = ephysin(i).stimulus(1,ikeep);
      ikeep2 = find(temp(:,1) >= ephysin(i).scanrange(1) & ...
                    temp(:,2) <= ephysin(i).scanrange(2));
      if (length(ikeep2) < size(temp,1))
          itoss = setdiff(1:size(temp,1),ikeep2);
          warnstr = sprintf(['File %s: file boundaries exceeded by [%g, %g] seconds;\n' ...
              'some valve transitions tossed because intervals exceed recording range'],...
              ephysin(i).basefilename, ...
              (ephysin(i).scanrange(1) - min(temp(:,1)))/ephysin(i).scanrate,...
              (max(temp(:,2))-ephysin(i).scanrange(2))/ephysin(i).scanrate);
          warning(warnstr);
      end
      intervals{i} = temp(ikeep2,:);
      identity{i} = identity{i}(ikeep2);
    else
      ikeep2i = find(ikeep < size(ephysin(i).stimulus,2));
      ikeep2 = ikeep(ikeep2i);
      temp = [ephysin(i).stimulus(2,ikeep2)',ephysin(i).stimulus(2,ikeep2+1)'];
      intervals{i} = temp;
      identity{i} = ephysin(i).stimulus(1,ikeep2);
    end
  end
%  if (nfiles == 1)
%    identity = identity{1};
%    intervals = intervals{1};
%  end
