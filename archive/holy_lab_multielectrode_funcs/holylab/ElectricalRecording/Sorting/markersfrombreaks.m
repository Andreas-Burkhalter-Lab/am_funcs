function [markers, labels] = markersfrombreaks(...
  cellSortFile,chanfile,files,data);

% Way to extract the information entered as markers
% for all cells
%
% "markers" contains the time in continuous time (ie, 
% time reference to the beginning of the first file, not
% the beginning of the file in which that file is found)
%
% not very carefully written; needs revision (esp
% a way to warn the user if they've chosen suspicious
% combinations of start and stop points)
%
% RCH 2006-01-20 wrote it

% get the information out of the file & calculate stuff we need
breakpointinfo = getbreakpoints(cellSortFile,chanfile);
if isstr(breakpointinfo)
  if strcmp(breakpointinfo,'none')
    intervals = 'none';
    return
  end
end
  
nFiles = length(files);
nBreaks = length(breakpointinfo);

% Make sure know correct offsets (in scans - note that as currently written
% this requires all files to have the same scanrate) for each filestart...
scanrate = data(1).header.scanrate;
for nthFile = 1:nFiles
  if data(nthFile).header.scanrate ~= scanrate
    error('Scanrate changes between files; this function not written to handle that')
  end
end
% create a list of offsets for each file so that at the end we can 
% turn our interval times into REAL across-file times, instead of the
% hack job (which ignores time gaps between files) we've been using...
d = data(1).header.date;
t = data(1).header.time;
realStartInDays = datenum([d t]);
for nthFile = 1:nFiles
  d = data(nthFile).header.date;
  t = data(nthFile).header.time;
  startInDays = datenum([d t]);
  offsetInDays = startInDays-realStartInDays;
  offsetInScans(nthFile) = offsetInDays*24*60*60*scanrate;
end

%create a list of adjusted scanrages so we can compare breaktimes
%across files
scanranges = NaNs(nFiles,2);
for nthFile = 1:nFiles
  scanranges(nthFile,:) = data(nthFile).scanrange;
end
adjustedscanranges = scanranges;
for nthFile = 2:nFiles
  adjustedscanranges(nthFile,:) = ...
    adjustedscanranges(nthFile,:) + offsetInScans(nthFile);
end
%add adjusted versions of timing to breakpointinfo
for nthBreak = 1:nBreaks
  fileNum = breakpointinfo(nthBreak).fileIndex;
  if fileNum == 1
    breakpointinfo(nthBreak).contTime = breakpointinfo(nthBreak).fileTime;
  else
    breakpointinfo(nthBreak).contTime = ...
      breakpointinfo(nthBreak).fileTime + ...
      offsetInScans(fileNum);
  end
end
%make a list of breaks indexed like breakpoint info that
%lets you know whether each breakpoint is a general marker or not
b = NaNs(1,nBreaks);
comments = cell(1,nBreaks);
times = NaNs(1,nBreaks);
for nthBreak = 1:nBreaks
  breakInfo = breakpointinfo(nthBreak);
  if strcmp(breakInfo.action,'mark_all') 
    b(nthBreak) = 1;
    comments{nthBreak} = breakInfo.comment;
    times(nthBreak) = breakInfo.contTime;
  else
    b(nthBreak) = 0;
    comments{nthBreak} = [];
    times(nthBreak) = breakInfo.contTime;
  end
end

ikill = find(b==0);
if length(ikill>0)
  nkill = length(ikill);
  times(ikill) = [];
  for nthKill = 1:nkill
    comments{ikill(nthKill)} = [];
  end
end

markers = times;
labels = comments;

% %in order to do the next step, we need to have these arranged
% %in order of time, so create a conversion matrix
% times = fieldout(breakpointinfo,'contTime');
% [y,iSorted] = sort(times);
%     % now y = times(iSorted)
% bSorted = b(iSorted);
% %edit the list to eliminate "redundant" breakpoints
% %(eg an all start and then a cell start - only a start
% %where the next breakpoint is a stop or visa versa 
% %will be used; endpoints are exceptions); also figure
% %out if 0 should be used as a start point or not and
% %whether the end point should be used as a stop
% rel = find( abs(bSorted)>0 );
% nRel = length(rel);
% if nRel == 0
%   intervals = 'none'
%   return
% end
% if bSorted(rel(1)) == -1
%   usestart = 1;
% else
%   usestart = 0;
%   if nRel > 1
%     if bSorted(rel(2)) == 1
%       b( find(iSorted == (rel(1)) ) ) = 0;
%     end
%   end
% end
% for nthRel = 2:(nRel-1)
%   this = bSorted(rel(nthRel));
%   prev = bSorted(rel(nthRel-1));
%   next = bSorted(rel(nthRel+1));
%   if this == 1
%     if next == 1
%       b( find(iSorted == (rel(nthRel)) ) ) = 0;
%     end
%   elseif this == -1
%     if prev == -1
%       b( find(iSorted == (rel(nthRel)) ) ) = 0;
%     end
%   else
%     error('something is off in indexing of relevant breakpoints')
%   end
% end
% if bSorted(rel(end)) == 1
%   usestop = 1;
% else
%   usestop = 0;
%   if nRel > 1
%     if bSorted(rel(end-1)) == -1
%       b( find(iSorted == (rel(end)) ) ) = 0;
%     end
%   end
% end
% %ok, now, for each file break we have to figure out if 
% %its start point is included in an "on period",
% %its stop point is included in an "on period", or neither
% nFilebreaks = nFiles-1;
% includeFbreak = NaNs(1,nFilebreaks);
% for nthFbreak = 1:nFilebreaks
%   fileNums = fieldout(breakpointinfo,'fileIndex');
%   ibefore = find((fileNums < (nthFbreak+1))&(abs(b)>0));
%   if ((length(ibefore)<1) & (usestart == 0))
%     includeFbreak(nthFbreak) = 0;
%   elseif ((length(ibefore)<1) & (usestart == 1))
%     includeFbreak(nthFbreak) = 1;
%   else
%     timesBefore = fieldout(breakpointinfo(ibefore),'contTime');
%     bBefore = b(ibefore);
%     [timesBefSorted,iSort] = sort(timesBefore);
%     bBefSorted = bBefore(iSort);
%     lastB = bBefSorted(end);
%     if lastB == -1
%       includeFbreak(nthFbreak) = 0;
%     elseif lastB == 1
%       includeFbreak(nthFbreak) = 1;
%     end
%   end
% end
% endpointcompilation = [usestart includeFbreak usestop];
% %finally, go through and compile a list, organized 
% %by file number, of start and stop times
% intervals = cell(1,nFiles);
% for nthFile = 1:nFiles
%   scanrange = data(nthFile).scanrange;
%   starts = [];
%   stops = [];
%   if endpointcompilation(nthFile) == 1
%     starts = [starts scanrange(1)];
%   end
%   if endpointcompilation(nthFile+1) == 1;
%     stops = [stops scanrange(2)];
%   end
%   fileNums = fieldout(breakpointinfo,'fileIndex');
%   iStarts = find( (fileNums==nthFile) & (b==1) );
%   startTimes = fieldout(breakpointinfo(iStarts),'fileTime');
%   starts = [starts startTimes];
%   iStops = find( (fileNums==nthFile) & (b==-1) );
%   stopTimes = fieldout(breakpointinfo(iStops),'fileTime');
%   stops = [stops stopTimes];
%   starts = sort(starts);
%   stops = sort(stops);
%   intervals{nthFile} = [starts' stops'];
% end
% intervals % for now want screen output to check
% 
