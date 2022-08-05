function h = msplot(scan,plotparams)
  if (nargin < 2)
    plotparams = struct;
  end
  deftype = 'time';
  if isfield(scan,'percentB')
    deftype = 'gradient';
  end
  plotparams = default(plotparams,...
    'type','image',...
    'xaxis',deftype, ...
    'tag_thresh',0.01,...
    'trange',[0 Inf],...
    'gradrange',[0 100],...
    'log',false,...
    'sqrt',false,...
    'show_xlabel',true);
  
  if ~isfield(plotparams,'scan_index')
    t = [scan.scan_time]/60;  % Convert to minutes
    sflag = (t >= plotparams.trange(1) & t <= plotparams.trange(2));
    if isfield(scan,'percentB')
      pB = [scan.percentB];
      sflag = sflag & (pB >= plotparams.gradrange(1) & pB <= plotparams.gradrange(2));
    end
  else
    % explicitly supply scan index
    sflag = false(1,length(scan));
    sflag(plotparams.scan_index) = true;
  end
  switch plotparams.type
    case 'mz'
      h = plot([scan(sflag).mz]);
      axis tight
      yl = get(gca,'YLim');
      yl = mean(yl) + [-1; 1]*diff(yl)/4;
      x = [0 cumsum([scan(sflag).n_points])];
      line([x; x], repmat(yl(:),1,length(x)),'Color','r')
      if plotparams.show_xlabel
        xlabel('Point index')
      end
      ylabel('m/z')

    case 'totI'
      sflag = sflag & ([scan.ms_level] == 1);
      t = t(sflag);
      I = [scan(sflag).totIntensity];
      h = plot(t,I);
      if plotparams.show_xlabel
        xlabel('Time (s)')
      end
      ylabel('Intensity')
   
    case 'image'
      sflag = sflag & ([scan.ms_level] == 1);
      interpflag = false;
      if strcmp(plotparams.xaxis,'time')
        % Interpolate values at evenly-spaced times
        scan_x = t(sflag);
        interpflag = true;
      elseif strcmp(plotparams.xaxis,'gradient')
        % Only include the first monotonic region of the gradient
        scan_x = [scan.percentB];
        negflag = [false, scan_x(2:end) < scan_x(1:end-1)];
        cnegflag = cumsum(double(negflag)) > 0;
        sflag = sflag & ~cnegflag;
        scan_x = [scan(sflag).percentB];
        interpflag = true;
      end
      [M,minmax] = msdecimate(scan(sflag),plotparams);
      
      % Convert scan # to times or gradient (and vice versa)
      n_scans = size(M,2);  % might be smaller than the original
      x_range = [1 n_scans];
      if interpflag
        x_evenly_spaced = linspace(scan_x(1),scan_x(end),n_scans);
        scan_spacing = interp1(scan_x,1:n_scans,x_evenly_spaced);
        M = interp1(M',scan_spacing)';
        x_range = scan_x([1 end]);
      end
      
      if plotparams.log
        M = log10(M+1);
      end
      if plotparams.sqrt
        M = sqrt(M);
      end

      h = imagesc(x_range,minmax,M);
      set(gca,'YDir','normal','XLim',x_range,'YLim',minmax)
      ylabel('m/z')
      if plotparams.show_xlabel
        if strcmp(plotparams.xaxis,'time')
          xlabel('Time (min)')
        elseif strcmp(plotparams.xaxis,'gradient')
          xlabel('% B')
        else
          xlabel('msLevel=1 scan index')
        end
      end
      set(gca,'TickDir','out')

    case 'slice-scan'
      mz = scan(plotparams.scan_index).mz;
      I = scan(plotparams.scan_index).intensity;
      if (scan(plotparams.scan_index).ms_level > 1)
        h = stem(mz,I,'.');
      else
        h = plot(mz,I);
      end
      if plotparams.tag_thresh < 1
        mx = max(I);
        ops = struct('polarity',1,'isMergeAdjacent',true,'adjacentNScans',length(I)/20);
        [pkIndex,pkVal] = findpeaks_multichan(I,plotparams.tag_thresh*mx,ops);
        for i = 1:length(pkIndex)
          thismz = double(mz(pkIndex(i)));
          text(thismz,double(pkVal(i)),sprintf('%0.2f',thismz));
        end
      end
      if plotparams.show_xlabel
        xlabel('m/z')
      end
      ylabel('Ion current')
      set(gca,'TickDir','out')

    case 'sum-scan'
      sflag = sflag;
      [mzI,scanIndex,I,mzLookup] = mscollect_points(scan(sflag),struct('mzI',true));
      Itot = accumarray(mzI(:),I(:));
      h = plot(mzLookup,Itot);
      if plotparams.show_xlabel
        xlabel('m/z')
      end
      ylabel('Total ion current')
      set(gca,'TickDir','out')
  end
end