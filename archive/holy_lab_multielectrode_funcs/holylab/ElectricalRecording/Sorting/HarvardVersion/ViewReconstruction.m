function ViewReconstruction(snips,tsnips,clustindx,sniprange,trange)
fview = figure('Position',[21   475   961   231],'Visible','off','BackingStore','off');        % Don't show until done
%fslide = figure('Position',[21 341 961 100]);
if (nargin < 5)
        trange = [0,tsnips(end)+sniprange(2)];
end
ncells = length(clustindx);
% Reconstruct waveforms
for i = 1:ncells
        [wcell{i},tcell{i}] = ReconstructWaveform(snips(:,clustindx{i}),tsnips(clustindx{i}),sniprange,trange);
end
xlim = [0 max(tsnips)+sniprange(2)];
ylim = [min(min(snips)) max(max(snips))];
set(gca,'NextPlot','add','XLim',xlim,'YLim',ylim);
co = get(gca,'ColorOrder');
ncol = size(co,1);
for i = 1:ncells
        plot(tcell{i},wcell{i},'Color',co(mod(i-1,ncol)+1,:),'HitTest','off','EraseMode','none');
end
ylabel('Voltage (V)');
xlabel('Time (scan #)');
sliderwindow(gca);
set(fview,'Visible','on');
