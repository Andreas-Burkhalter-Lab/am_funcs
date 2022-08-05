function imsaveas(h,filename,format)
  % Work around matlab bug in SAVEAS
  imhandle = findobj(h,'type','image');
  im = get(imhandle,'CData');
  clim = get(get(imhandle,'Parent'),'CLim');
  slope = 1/diff(clim);
  imrs = slope*(im-clim(1));
  imwrite(imrs,filename,format)
  
    