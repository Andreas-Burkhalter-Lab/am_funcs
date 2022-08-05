function ephysout = ephyssubrange(ephysin,intervals)
% EPHYSSUBRANGE: chop ephys data into time intervals
%  ephysout = ephyssubrange(ephysin,intervals)
%
% where
%   ephysin is the input ephys data (a structure, see EPHYS)
%   intervals specifies to time intervals to use in chopping up the data.
%     It is an nintervals-by-2 matrix, with column 1 giving start times
%     and column 2 giving end times, in scans.
%
% ephysout is a structure array of ephys data.
%
% Note also that ephysin may be a structure array, and intervals a cell
%   array of interval matrices.  In such cases, ephysout is a cell array
%   of ephys structure arrays organized in the same way as ephysin.  In
%   most cases you'll want to concatenate these into a single structure
%   array, using cat(2,ephysout{:}), but in some cases you may want to
%   preserve the input structure.
%
% Check to see if the sorted-cell stuff really works!
%
% See also: EPHYS, EPHYSFETCH, EPHYSSUBCHAN, EPHYSSUBCELL, EPHYSTAG, INTERVALSFROMSTIM.
  
% Changelog:
%   TEH, 05-03-2004: generalize to work with arrays
%   RCH, 11-09-2005: added ability to handle files with no intervals
%   requested
  
  % Supply default values
  if (nargin < 2)
    ephysout = ephysin;
    return;
  end
  
  if ~iscell(intervals)
    intervals = {intervals};
  end
  
  nepochsi = length(ephysin(:));
  if (length(intervals(:)) ~= nepochsi)
    error(['The number of interval matrices does not agree with the number ' ...
           'of elements in ephysin.']);
  end
  % Loop over the elements in ephysin
  ephysout = cell(size(ephysin));
  for i = 1:nepochsi
    if isempty(intervals{i})
      ephysout{i} = [];
    else
      ephysout{i} = ephyssubrangework(ephysin(i),intervals{i});
    end
  end
  
  if (nepochsi == 1)
    ephysout = ephysout{1};
  end
  
  
% This next function does the real work. It was written first, and the
% code to generalize to arrays was written later, and it was easier to
% make this a subfunction.
function ephysout = ephyssubrangework(ephysin,intervals)  
  % Check that the given intervals are within the range contained
  % in the input data
  if (min(intervals(:,1)) < ephysin.scanrange(1) | ...
      max(intervals(:,2)) > ephysin.scanrange(2))
    error('Desired intervals fall outside of range of data');
  end

  % Check that we aren't trying to get a subinterval on
  % snippet or spikepeak data without having spiketimes available
  if ((isfield(ephysin,'snippets') | isfield(ephysin,'snippeaks') | ...
       isfield(ephysin,'snipindex')) & ~isfield(ephysin,'sniptimes'))
    error('Trying to get a subinterval with no spike timing data available')
  end


  
  % We can get started!
  nintervals = size(intervals,1);
  nchannels = length(ephysin.channels);
  ncells = 0;
  if isfield(ephysin,'celltimes')
    ncells = length(ephysin.celltimes);
  end
  % Copy over any unmodified data
  fieldstomodify = {'stimulus','wave','envelope','snippets', ...
                    'sniptimes','snippeaks','snipindex','celltimes', ...
                    'scanrange','cellscanrange'};
  ephysbase = copyotherfields(ephysin,fieldstomodify);
  % Pre-allocate storage for cell arrays
  snipnames = {'sniptimes','snipindex','snippets','snippeaks'};
  for i = 1:length(snipnames)
    if isfield(ephysin,snipnames{i})
      ephysbase.(snipnames{i}) = cell(1,nchannels);
    end
  end
  cellnames = {'celltimes','cellscanrange'};
  for i = 1:length(cellnames)
    if isfield(ephysin,cellnames{i})
      ephysbase.(cellnames{i}) = cell(1,ncells);
    end
  end
  ephysout = repmat(ephysbase,1,nintervals);
  
  %
  % Subrange the modified fields
  %
  % scanrange
  for i = 1:nintervals
    ephysout(i).scanrange = intervals(i,:);
  end

  % stimulus 
  if isfield(ephysin,'stimulus')
    for i = 1:nintervals
      lindx = find(ephysin.stimulus(2,:) <= intervals(i,1));
      rindx = find(ephysin.stimulus(2,:) <= intervals(i,2));
      rng = [lindx(end):rindx(end),rindx(end)];
      ephysout(i).stimulus = ephysin.stimulus(:,rng);
      ephysout(i).stimulus(2,[1 end]) = intervals(i,:);
    end
  end
  
  % wave
  if isfield(ephysin,'wave')
    intervals_shift = intervals - ephysin.scanrange(1) + 1;
    for i = 1:nintervals
      ephysout(i).wave = ephysin.wave(:, ...
        intervals_shift(i,1):intervals_shift(i,2));
    end
  end

  % envelope
  if isfield(ephysin,'envelope')
    intervals_shift = round((intervals - ephysin.scanrange(1)) / ...
              ephysin.envelopedecimate) + 1; 
    for i = 1:nintervals
        ephysout(i).envelope = ephysin.envelope(:, ...
          intervals_shift(i,1):intervals_shift(i,2));
    end
  end
  
  % sniptimes, snipindex, snippets, snippeaks are relative to the times in
  % sniptimes, so do all these together
  if isfield(ephysin,'sniptimes')
    for j = 1:nchannels
      keepindex = timesininterval(ephysin.sniptimes{j},intervals);
      for k = 1:nintervals
        ephysout(k).sniptimes{j} = ephysin.sniptimes{j}(keepindex{k});
        if isfield(ephysin,'snipindex')
          ephysout(k).snipindex{j} = ephysin.snipindex{j}(keepindex{k});
        end
        if isfield(ephysin,'snippets')
          ephysout(k).snippets{j} = ephysin.snippets{j}(:,keepindex{k});
        end
        if isfield(ephysin,'snippeaks')
          ephysout(k).snippeaks{j} = ephysin.snippeaks{j}(keepindex{k});
        end
      end  % over intervals
    end % over channels
  end % if isfield(ephysin,'sniptimes')
  
  % cells
  if isfield(ephysin,'celltimes')
    for j = 1:ncells
      if isfield(ephysin,'cellscanrange')
        this_cell_intervals = ephysin.cellscanrange{j};
      else
        this_cell_intervals = ephysin.scanrange;
      end
      for cellIntervalIndex = 1:size(this_cell_intervals,1)
        [cell_intervals,intervalRetained] = IntersectIntervals(...
          this_cell_intervals(cellIntervalIndex,:),intervals);
        if isempty(cell_intervals)
          continue;
        end
        keepindex = timesininterval(ephysin.celltimes{j},cell_intervals);
        intervalRetained = find(intervalRetained);
        for k = 1:length(intervalRetained)
          indxR = intervalRetained(k);
          if isempty(ephysout(indxR).celltimes{j})
            ephysout(indxR).celltimes{j} = ephysin.celltimes{j}(keepindex{k});
          else
            ephysout(indxR).celltimes{j} = [ephysout(indxR).celltimes{j}; ...
              ephysin.celltimes{j}(keepindex{k})];
          end
        end
      end
    end
  end
  
  if isfield(ephysin,'cellscanrange')
    for k = 1:nintervals
      for j = 1:ncells
        ephysout(k).cellscanrange{j} = IntersectIntervals(intervals(k,:),...
          ephysin.cellscanrange{j});
      end
    end
  end
  
  return
  
  %if isfield(ephysin,'cells')
  %  intervals_units = intervals_secs/ephysin.cells.tosecs;
  %  ephysout(i).cells = copyotherfields(ephysin.cells,'data');
  %  if isfield(ephysin.cells,'data')
  %    ephysout(i).cells.data = cell(1,length(ephysin.cells.data));
  %  end
  %end
  %    if isfield(ephysin.spiketimes,'data')
  %      for j = 1:nchannels
  %        keepindex = timesininterval(ephysin.spiketimes.data{chanIndex(j)}, ...
  %          intervals_units);
%           for k = 1:nintervals
%             ephysout(k).spiketimes.data{j} = ...
%               ephysin.spiketimes.data{chanIndex(j)}(keepindex{k});




