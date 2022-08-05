function hax = plot_stack_row(im)
% PLOT_STACK_ROW: plot images in a stack along a row
% Syntax:
%   hax = plot_stack_row(im)
% where
%   im is a m-by-n-by-K array of K images
% and
%   hax is a vector of handles to the output axes.
  
  nplots = size(im,3);
  figure
  hax = zeros(1,nplots);
  for i = 1:nplots
    hax(i) = subplot(1,nplots,i);
    imshowsc(im(:,:,i))
  end
  clim = [min(im(:)) max(im(:))];
  set(hax,'CLim',clim,'XTick',[],'YTick',[]);
end
