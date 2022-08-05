function intervals = timemarkers2intervals(timeMarker,nscans,n_cells)
% TIMEMARKERS2INTERVALS: parse breakpoints from CASS
% Syntax:
%   intervals = timemarkers2intervals(timeMarker,nscans,n_cells)
% where
%   timeMarker is a structure array of breakpoints, of the type saved in
%     the sort_info structure output by CASS.
%   nscans is a vector containing the total number of scans in a set of
%     files (e.g., from [sortheader.nscans])
%   n_cells is the total number of cells on this CASS channel (can be
%     obtained from max(landmarkCluster)
% and
%   intervals is a n_cells-by-n_files cell array, each element of which
%     is an n_intervals-by-2 matrix giving the range of valid scan #s for
%     the given cell,file combination.
%
% See also: INTERVALSFROMBREAKS.
  
% Basically, this differs from INTERVALSFROMBREAKS only in not being so
% file-dependent; it processes from information in memory.

% Copyright 2007 by Timothy E. Holy
  
  n_files = length(nscans);
  intervals = cell(n_cells,n_files);
  % Make a quick exit if no timeMarkers
  if isempty(timeMarker)
    % Supply default intervals: full file
    for cellIndex = 1:n_cells
      for fileIndex = 1:n_files
        intervals{cellIndex,fileIndex} = [1 nscans(fileIndex)];
      end
    end
    return;
  end
  % Parcel out relevant time markers by cell #. Use a matrix:
  %      [globaltime (cumulative scan #)
  %       eventtype (-1 for stop, 1 for start)]
  scan_offset = [0 cumsum(nscans)];
  cellevents = cell(1,n_cells);
  for cellIndex = 1:n_cells
    cellevents{cellIndex} = zeros(0,2);
  end
  for eventIndex = 1:length(timeMarker)
    tm = timeMarker(eventIndex);
    % Limit event times to within files
    fileTime = min(tm.fileTime,nscans(tm.fileIndex));
    globalTime = fileTime + scan_offset(tm.fileIndex);
    switch(tm.action)
      case 'start_cell',
        if tm.clustNum > 0
          cellevents{tm.clustNum}(end+1,:) = [globalTime 1];
        end
      case 'stop_cell',
        if tm.clustNum > 0
          cellevents{tm.clustNum}(end+1,:) = [globalTime -1];
        end
      case 'start_all',
        for cellIndex = 1:n_cells
          cellevents{cellIndex}(end+1,:) = [globalTime 1];
        end
      case 'stop_all',
        for cellIndex = 1:n_cells
          cellevents{cellIndex}(end+1,:) = [globalTime -1];
        end
      otherwise
        continue %irrelevant marker type
    end
  end
  % Put events in order
  for cellIndex = 1:n_cells
    cellevents{cellIndex} = sortrows(cellevents{cellIndex});
  end

  % Now the final processing: loop over the individual cells, polish up the
  % set of events (put in sentinels, remove redundancies), and break into
  % intervals
  fake_events = [1 1; scan_offset(end) -1];  %start/stop at extremes
  for cellIndex = 1:n_cells
    % Put in sentinel points to make it easier to process start/stop info
    cur_cellevents = cellevents{cellIndex};
    if ~isempty(cur_cellevents)
      if (cur_cellevents(1,2) == -1)
        % First event was a stop, user must want to include beginning
        cur_cellevents = [fake_events(1,:); cur_cellevents];
      end
      if (cur_cellevents(end,2) == 1)
        % Last event was a start, user must want to continue until end
        cur_cellevents = [cur_cellevents; fake_events(2,:)];
      end
    else
      % No events detected, insert start/stop at beginning/end of recording
      cur_cellevents = fake_events;
    end
    % Go through and remove redundancies, e.g., a stop followed by
    % another stop. We always choose to keep the first of these redundant
    % events
    killIndex = find(cur_cellevents(2:end,2) == cur_cellevents(1:end-1,2));
    cur_cellevents(killIndex+1,:) = [];
    % OK, now we can count on the fact that events are ordered in the
    % following way: start, stop, start, stop, etc.
    cur_intervals = reshape(cur_cellevents(:,1),...
      [2 size(cur_cellevents,1)/2])';
    % Intersect these with file scan ranges
    for fileIndex = 1:n_files
      cur_scans = scan_offset([0 1]+fileIndex) + [1 0];
      intersect_int = IntersectIntervals(cur_scans,cur_intervals);
      % Go back to real scan units, i.e. remove the scan offset
      intervals{cellIndex,fileIndex} = intersect_int - scan_offset(fileIndex);
    end
  end
