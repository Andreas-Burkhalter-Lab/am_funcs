function msbrowse(scan)
% MSBROWSE: interactively inspect a mass-spec file
%
% Syntax:
%   msbrowse(scan)
% where
%   scan is a structure array represent a mass-spec data file, for example
%     of the type returned by msload_mzXML.
%
% This plots the intensity as a function of elution time (along the x-axis)
% and m/z (along with y-axis). When browsing,
%   left-click-drag selects a zoom region
%   middle-click returns to the full view
%   right-click plots the intensity at a particular time point, as well as
%     any MS/MS data available at that time (in two other figure windows)
%   right-click-drag selects a region of time over which to average, and
%     plots the intensity
%
% See also: MSLOAD_MZXML.

% Copyright 2009-2010 by Timothy E. Holy

  hfprof = figure('Name','profile');
  hfmsms = figure('Name','MS/MS');
  hfig = figure;
  himg = msplot(scan,struct('type','image'));
  hax = gca;
  set(hfig,'Units','normalized');
  setappdata(hfig,'scan',scan);
  setappdata(hfig,'hfprof',hfprof);
  setappdata(hfig,'hfmsms',hfmsms);
  set(hax,'NextPlot','replacechildren','ButtonDownFcn',@msbclick);
  set(himg,'HitTest','off')
end

function msbclick(src,eventdata)
  hax = src; %get(src,'Parent');
  hfig = get_parent_fig(hax);
  cp = get(hax,'CurrentPoint');
  selType = get(hfig,'SelectionType');
  scan = getappdata(hfig,'scan');
  hfprof = getappdata(hfig,'hfprof');
  hfmsms = getappdata(hfig,'hfmsms');
  switch lower(selType)
    case 'normal'
      [x,y] = msb_get_drag_region(hax);
      if (range(x) == 0 || range(y) == 0)
        return
      end
      ops = struct('type','image',...
        'min_mz',y(1),...
        'max_mz',y(2),...
        'trange',x);
      axes(hax);
      himg = msplot(scan,ops);
      set(himg,'HitTest','off')
    case 'alt'
      t = msb_get_drag_region(hax);
      flag = [scan.ms_level] == 1;
      scanMap = find(flag);  % to convert back to all scans (not just ms_level==1)
      tscan = [scan(flag).scan_time]/60;
      tIndex = find(tscan >= t(1) & tscan <= t(2));  % the selected scans
      if (length(tIndex) <= 1)
        % Plot the profile and MS/MS
        % In case tIndex is empty, let's just find the closest one
        [tmp,tIndex] = min(abs(t(1)-tscan));
        % Convert the time to a scan index
        scanIndex = scanMap(tIndex);
        % Plot the profile data
        figure(hfprof)
        h = msplot(scan,struct('type','slice-scan','scan_index',scanIndex));
        title(sprintf('Time: %g',t(1)))
        set(h,'HitTest','off');
        set(gca,'ButtonDownFcn',@msbintegrate);
        % Plot the MS/MS data
        figure(hfmsms)
        indx = scanMap(tIndex)+1:scanMap(tIndex+1)-1;
%         indx_base = -5:5;
%         indx = scanIndex + indx_base;
%         ms_level = [scan(indx).ms_level];
%         indx = indx(ms_level > 1);
        precursor_mz = [scan(indx).precursor_mz];
        [umz,tmp,pcIndex] = unique(precursor_mz);
        clf
        nmsms = length(umz);
        xlimm = zeros(nmsms,2);
        msmsax = zeros(1,nmsms);
        for i = 1:nmsms
          msmsax(i) = subplot(nmsms,1,i);
%           thismz = find(pcIndex == i);
%           [mnoffset,closest] = min(abs(indx_base(thismz)));
%           msmsIndex = indx(thismz(closest));
          msmsIndex = indx(pcIndex == i);
          msplot(scan,struct('type','slice-scan','scan_index',msmsIndex,'show_xlabel',i == nmsms));
          title(sprintf('Time: %g, Precurs %g, Precurs intensity %g',t(1),umz(i),scan(msmsIndex).precursor_I))
          xlimm(i,:) = get(gca,'XLim');
        end
        set(msmsax,'XLim',max(xlimm,[],1));
        set(msmsax(1:end-1),'XTick',[]);
      else
        % Average over a temporal region and plot the intensity
        figure(hfprof)
        msplot(scan,struct('type','sum-scan','scan_index',scanMap(tIndex)));
      end
    case 'extend'
      ops = struct('type','image');
      axes(hax)
      himg = msplot(scan,ops);
      set(himg,'HitTest','off')
  end
end

function msbintegrate(src,eventdata)
  hax = src; %get(src,'Parent');
  hfig = get_parent_fig(hax);
  cp = get(hax,'CurrentPoint');
  selType = get(hfig,'SelectionType');
  switch selType
    case 'alt'
      [pos,linehandle,vertexIndex] = findpoint(cp(1,:),hax);
      ydata = get(linehandle,'YData');
      rng = max(1,vertexIndex-5):min(length(ydata),vertexIndex+5);
      fprintf('m/z: %0.4f, integral: %g\n',pos(1),sum(ydata(rng)));
  end
end

function [x,y] = msb_get_drag_region(hax)
  % Start dragging to zoom on region
  finalRect = rbbox;
  xn = [0 finalRect(3)] + finalRect(1);
  yn = [0 finalRect(4)] + finalRect(2);
  [x,y] = norm2data(xn,yn,hax);
end
