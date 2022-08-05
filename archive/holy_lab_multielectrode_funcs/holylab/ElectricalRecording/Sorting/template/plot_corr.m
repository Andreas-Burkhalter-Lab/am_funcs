function plot_corr(hAxes, spiketimes, options)
% plot_corr(hAxes, spiketimes, options)

    options = default(options,'trange', 5000); % unit is 0.1ms if scanrate is 10k
    
    trange=options.trange;
    tpos=getpixelposition(hAxes);
    npix=tpos(3);
    %nbins = ceil(npix/2);
    nbins = npix;
    
    isCrossCorr=iscell(spiketimes);
    
    if(isCrossCorr)
       n = crosscorrspike(spiketimes{1}, spiketimes{2}, trange, nbins);
       binwidth = 2*trange/nbins;
       x = linspace(-trange+binwidth/2,trange-binwidth/2,nbins);
    else
       n = autocorrspike(spiketimes,trange,nbins);
       binwidth = trange/nbins;
       x = linspace(binwidth/2,trange-binwidth/2,nbins);
    end

    hbar = bar(x,n,1);
    if isfield(options,'col')
      set(hbar,'FaceColor',options.col,'EdgeColor','none');
    else
      set(hbar,'FaceColor','k');
    end

    if(isCrossCorr)
       set(hAxes,'XLim',[-trange trange])
    else
       set(hAxes,'XLim',[0 trange])
    end

    