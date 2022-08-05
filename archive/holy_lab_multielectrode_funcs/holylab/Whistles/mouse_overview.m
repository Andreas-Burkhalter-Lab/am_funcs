function [fls,nwhis,tsnif] = mouse_overview(mname)
% MOUSE_OVERVIEW: a function to summarize a mouse's response to diff. stim.
% Syntax:
%   mouse_overview(mousename)
%   [fls,nwhis,tsnif] = mouse_overview(mousename);
% where
%   mousename is a string containing the name assigned to a given mouse;
%     it is assumed that this name is at the beginning of each filename
%     for data collected on this mouse (e.g., m14fmu.bin would refer to a
%     mouse named 'm14', and be a test using female mouse urine (fmu).)
%
% The first syntax gives graphical output in the form of a bar chart. The
% second returns data on the files analyzed, the number of chirps in each
% file, and the total time spent sniffing.
  
% Copyright 2007 by Timothy E. Holy
  
  mintime = 0.015;  % Discard "snifs" lasting too little time (artifacts
                  % from putting the detector into the chamber and
                  % occasional noise-triggered events)
  duration = 209.76;
  snifwindow = 25;  % Window used in evaluating "sniff" fraction
  % Directory parsing
  dcur = pwd;
  sepIndex = strfind(mname,filesep);
  if ~isempty(sepIndex)
      dnew = mname(1:sepIndex(end));
      cd(dnew);
      mname = mname(sepIndex(end)+1:end);
  end
  fls = dirbyname([mname '*.bin']);
  flkeep = true(1,length(fls));
  for i = 1:length(fls)
      flkeep(i) = ~isempty(regexp(fls{i},[mname '\D']));
  end
  fls = fls(flkeep);
  fls = whistles_temporal_order(fls);
  for i = 1:length(fls)
    wt{i} = load([fls{i} '.whistime']);
    nwhis(i) = size(wt{i},1);
    [detdata,iv,h] = ReadSensor([fls{i} '.det'],...
                                struct('close','lo','fixend',250));
%     [detdata,iv,h] = ReadSensor([fls{i} '.det'],...
%                                 struct('close','hi','fixend',250));
    if (iv < 2.5)
        detdata = detdata(:,2:end);
    end
    dt = diff(detdata);
    indxsnif = find(dt > mintime);
    detdata = detdata(:,indxsnif);
    tsnif(i) = 0;
    if ~isempty(detdata)
      tstart = detdata(1,1);
      dt = detdata - tstart;
      indxl = find(dt(1,:) < snifwindow);
      tsnif(i) = sum(diff(detdata(:,indxl)))/snifwindow;
    end
    fprintf('%s %d %g\n',fls{i},nwhis(i),tsnif(i));
  end
  barh(tsnif*100)
  set(gca,'YTick',1:length(tsnif),'YTickLabel',fls,'YDir','reverse','XLim',[0 100],...
    'YLim',[0 length(tsnif)+1]);
  xlabel('% of time spent sniffing')
  cd(dcur)
  