function ClustSWCallback(action,hfig)
if (nargin == 1)
        hfig = gcbf;
end
switch(action)
case 'AutoCorr'
        selvec = getappdata(hfig,'selectflag');
        selclust = find(selvec);
        polygons = getappdata(hfig,'polygons');
        x = getappdata(hfig,'x');
        y = getappdata(hfig,'y');
        t = getappdata(hfig,'t');
        f = getappdata(hfig,'f');
        membership = ComputeMembership(x,y,polygons(selclust));
        indx = find(membership);
        for i = 1:length(t)
                findx = find(f(1,indx) == i);        % Find the subset in a given file
                tsub{i} = t{i}(f(2,indx(findx)));        % Pull out the times for this file
        end
        if (length(indx) > 0)
                autocorrfig(tsub,0.1,'s');
        else
                errordlg('Must select one or more clusters first');
        end
case 'Revert'
        %hax = findobj(hfig,'Tag','ClustAx');
        %axes(hax);
        %cla
        %set(gca,'Tag','ClustAx');
        setappdata(hfig,'polygons',getappdata(hfig,'polygons0'));
        clustnums = getappdata(hfig,'clustnums0');
        setappdata(hfig,'clustnums',clustnums);
        %setappdata(hfig,'hlines',[]);
        setappdata(hfig,'selectflag',zeros(size(clustnums)));
        if (strcmp(getappdata(hfig,'mode'),'density'))
                ClusterFunctions('DensityPlot',hfig);
        else
                ClusterFunctions('ScatterPlot',hfig);
        end
case 'Waveforms'
        wvfig = findobj(0,'Tag','WaveformSelectFig');
        if (isempty(wvfig) | ~ishandle(wvfig))
                wvfig = figure;
        end
        %set(wvfig,'KeyPressFcn','');
        cla
        x = getappdata(hfig,'x');
        y = getappdata(hfig,'y');
        wv = num2cell(getappdata(hfig,'wfms'),1);
        hlines = line([x;x],[y;y],'ButtonDownFcn','ClustSWCallback ShowWaveform','Color','r');
        set(hlines,{'UserData'},wv');
        setappdata(wvfig,'swylim',getappdata(hfig,'swylim'));
case 'ShowWaveform'
        h = 100; w = 180;
        ploc = get(0,'PointerLocation');
        figure('Position',[ploc(1),ploc(2)-h-10,w,h],'Resize','off');
        udata = get(gcbo,'UserData');
        plot(udata);
        set(gca,'YLim',getappdata(gcbf,'swylim'),'XLim',[1 length(udata)]);
otherwise
        error(['Do not recognize action ',action]);
end
