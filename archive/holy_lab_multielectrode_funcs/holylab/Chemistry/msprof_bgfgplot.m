function msprof_bgfgplot(hax,msprofdata)
  if ~ishandle(hax)
    msprofdata = hax;
    hax = gca;
  end
  % Draw the "error bars"
  hpbg = fillmm2((msprofdata.bg-msprofdata.bg_sem)',...
		 (msprofdata.bg+msprofdata.bg_sem)',...
		 msprofdata.mz','Parent',hax);
  set(hpbg,'FaceColor',[1 1 1]*0.6,'EdgeColor','none');
  hpfg = fillmm2((msprofdata.fg-msprofdata.fg_sem)',...
		 (msprofdata.fg+msprofdata.fg_sem)',...
		 msprofdata.mz','Parent',hax,'EdgeColor','none');
  set(hpfg,'FaceColor',[1 0.6 0.6]);
  % Draw the mean values
  hline = line(msprofdata.mz,[msprofdata.bg msprofdata.fg],...
	       'Parent',hax);
  set(hline,{'Color'},{[1 1 1]*0.2; [1 0 0]});
  
  set(hax,'TickDir','out') % so ticks don't get confused with peaks
end
  