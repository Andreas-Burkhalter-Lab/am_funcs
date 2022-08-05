function [intervals] = intervalsfrombreaks(...
  cellNumber,cellSortFile,chanfile,files,data,op);

% Way to turn breakpoint information into a set of 
% intervals relevant for a particular cell
%
% not very carefully written; needs revision (esp
% a way to warn the user if they've chosen suspicious
% combinations of start and stop points)
%
% RCH 2005-10-25 wrote it
% RCH 2006-02-13 fixed bug in handling of redundant breakpoints (doh!)
% RCH 2007-04-18 added ability to handle multiunit data
% RCH 2007-07-03 added check in case breakpoint in between files
% RCH 2007-07-05 MERGED TWO FORKED COPPIES! ACK!

if nargin < 6
    op = struct;
end
op = default(op,'show_intervals_as_running',0);
op = default(op,'move_breaks_past_fileend_to_fileend',1);

% get the information out of the file & calculate stuff we need
if iscell(cellSortFile)
    % means call is from something based around multiunit data, and it's
    % already figured out the appropriate breakpoint info for this cell for
    % us
    breakpointinfo = cellSortFile{cellNumber};
    if length(breakpointinfo) == 0
        breakpointinfo = 'none';
    end
else
    breakpointinfo = getbreakpoints(cellSortFile,chanfile);
end
if isstr(breakpointinfo)
  if strcmp(breakpointinfo,'none')
    intervals = 'none';
    return
  end
end
  
nFiles = length(files);
nBreaks = length(breakpointinfo);
%create a list of adjusted scanrages so we can compare breaktimes
%across files
scanranges = NaNs(nFiles,2);
for nthFile = 1:nFiles
  scanranges(nthFile,:) = data(nthFile).scanrange;
end
adjustedscanranges = scanranges;
for nthFile = 2:nFiles
  adjustedscanranges(nthFile,:) = ...
    adjustedscanranges(nthFile,:) + adjustedscanranges(nthFile-1,2);
end
%add adjusted versions of timing to breakpointinfo
for nthBreak = 1:nBreaks
  fileNum = breakpointinfo(nthBreak).fileIndex;
  if fileNum == 1
    breakpointinfo(nthBreak).contTime = breakpointinfo(nthBreak).fileTime;
  else
    breakpointinfo(nthBreak).contTime = ...
      breakpointinfo(nthBreak).fileTime + ...
      adjustedscanranges(fileNum-1,2);
  end
end
%make a list of breaks indexed like breakpoint info that
%lets you know whether for this cell, a given break is 
%an on, an off, or irrelveant
b = NaNs(1,nBreaks);
for nthBreak = 1:nBreaks
  breakInfo = breakpointinfo(nthBreak);
  if ( strcmp(breakInfo.action,'start_all') | ( strcmp(breakInfo.action,'start_cell') & breakInfo.clustNum == cellNumber))
    b(nthBreak) = 1;
  elseif ( strcmp(breakInfo.action,'stop_all') | ( strcmp(breakInfo.action,'stop_cell') & breakInfo.clustNum == cellNumber))
    b(nthBreak) = -1;
  else
    b(nthBreak) = 0;
  end
end
%in order to do the next step, we need to have these arranged
%in order of time, so create a conversion matrix
times = fieldout(breakpointinfo,'contTime');
[y,iSorted] = sort(times);
    % now y = times(iSorted)
bSorted = b(iSorted);
%edit the list to eliminate "redundant" breakpoints
%(eg an all start and then a cell start - only a start
%where the next breakpoint is a stop or visa versa 
%will be used; endpoints are exceptions); also figure
%out if 0 should be used as a start point or not and
%whether the end point should be used as a stop
rel = find( abs(bSorted)>0 );
nRel = length(rel);
if nRel == 0
  intervals = 'none'
  return
end
if bSorted(rel(1)) == -1
  usestart = 1;
else
  usestart = 0;
  if nRel > 1
    if bSorted(rel(2)) == 1
      b(iSorted(rel(1))) = 0;
    end
  end
end
for nthRel = 2:(nRel-1)
  this = bSorted(rel(nthRel));
  prev = bSorted(rel(nthRel-1));
  next = bSorted(rel(nthRel+1));
  if this == 1
    if next == 1
      b(iSorted(rel(nthRel))) = 0;
    end
  elseif this == -1
    if prev == -1
      b(iSorted(rel(nthRel))) = 0;
    end
  else
    error('something is off in indexing of relevant breakpoints')
  end
end
if bSorted(rel(end)) == 1
  usestop = 1;
else
  usestop = 0;
  if nRel > 1
    if bSorted(rel(end-1)) == -1
      b(iSorted(rel(end))) = 0;
    end
  end
end
%ok, now, for each file break we have to figure out if 
%its start point is included in an "on period",
%its stop point is included in an "on period", or neither
nFilebreaks = nFiles-1;
includeFbreak = NaNs(1,nFilebreaks);
for nthFbreak = 1:nFilebreaks
  fileNums = fieldout(breakpointinfo,'fileIndex');
  ibefore = find((fileNums < (nthFbreak+1))&(abs(b)>0));
  if ((length(ibefore)<1) & (usestart == 0))
    includeFbreak(nthFbreak) = 0;
  elseif ((length(ibefore)<1) & (usestart == 1))
    includeFbreak(nthFbreak) = 1;
  else
    timesBefore = fieldout(breakpointinfo(ibefore),'contTime');
    bBefore = b(ibefore);
    [timesBefSorted,iSort] = sort(timesBefore);
    bBefSorted = bBefore(iSort);
    lastB = bBefSorted(end);
    if lastB == -1
      includeFbreak(nthFbreak) = 0;
    elseif lastB == 1
      includeFbreak(nthFbreak) = 1;
    end
  end
end
endpointcompilation = [usestart includeFbreak usestop];
%finally, go through and compile a list, organized 
%by file number, of start and stop times
intervals = cell(1,nFiles);
for nthFile = 1:nFiles
  scanrange = data(nthFile).scanrange;
  starts = [];
  stops = [];
  if endpointcompilation(nthFile) == 1
    starts = [starts scanrange(1)];
  end
  if endpointcompilation(nthFile+1) == 1;
    stops = [stops scanrange(2)];
  end
  fileNums = fieldout(breakpointinfo,'fileIndex');
  iStarts = find( (fileNums==nthFile) & (b==1) );
  startTimes = fieldout(breakpointinfo(iStarts),'fileTime');
  starts = [starts startTimes];
  iStops = find( (fileNums==nthFile) & (b==-1) );
  stopTimes = fieldout(breakpointinfo(iStops),'fileTime');
  stops = [stops stopTimes];
  starts = sort(starts);
  stops = sort(stops);
  intervals{nthFile} = [starts' stops'];
  
  if op.move_breaks_past_fileend_to_fileend
      % very simple dumb check... might need upgrading if complexities
      % arise
      scanrange = data(nthFile).scanrange;
      iTooLate = find(intervals{nthFile}>scanrange(2));
      intervals{nthFile}(iTooLate) = scanrange(2);
  end
  
end

if op.show_intervals_as_running
    intervals % for now want screen output to check
end

