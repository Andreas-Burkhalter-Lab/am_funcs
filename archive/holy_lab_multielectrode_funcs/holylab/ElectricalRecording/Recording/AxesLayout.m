function [ax,width] = AxesLayout (hfig,grid,axesType)

l = grid(1);
b = grid(2);
w = grid(3);
h = grid(4);

if (axesType ==1)
    load MEAInitialization.mat;
elseif(axesType == 2)
    load BundleInitialization.mat;
else
    load FourChannelInitialization.mat;
end

prevAxes = findobj(hfig,'Tag','ExperimentAxes');
prevAxes = [prevAxes findobj(hfig,'Tag','HiddenAxes')];
delete(prevAxes);


for j=1:prod(size(ax))
        left = l + ax(j).position(1) * w;
        bottom = b + ax(j).position(2) * h;
        width = ax(j).position(3) * w;
        height = ax(j).position(4) * h;
        
        h1 = axes('Parent',hfig, ...
                'Units','pixels', ...
                'CameraUpVector',[0 1 0], ...
                'Color',[1 1 1], ...
                'Position',[left,bottom,width,height], ...
                'Tag','ExperimentAxes', ...
        'Box','on', ...
        'ButtonDownFcn','RecCB ChannelStatus',...
        'NextPlot','replacechildren',...
        'XTick',[ ], ...
        'YTick',[ ], ...
                'XColor',[0 0 0], ...
                'YColor',[0 0 0], ...
                'ZColor',[0 0 0]);
        
        h2 = axes('Parent',hfig, ...
                'Units','pixels', ...
                'CameraUpVector',[0 1 0], ...
                'Color',[1 1 1], ...
                'Position',grid, ...
        'Visible','off',...
        'HandleVisibility','on',...
        'XTick',[ ], ...
        'YTick',[ ], ...
                'Tag','HiddenAxes', ...
                'XColor',[0 0 0], ...
                'YColor',[0 0 0], ...
                'ZColor',[0 0 0]);
    
    hAxes = findobj(hfig,'Tag','ExperimentAxes');

    ax(j).axhandle = hAxes(1,1);

        if((prod(size(ax))) == 64)
        ax(j).labelhandle = text(left-18,bottom+58,ax(j).label,'Parent',h2,'Units','pixels','Color','b');
        else
        ax(j).labelhandle = text(left-18,bottom+128,ax(j).label,'Parent',h2,'Units','pixels','Color','b');
        end
end

