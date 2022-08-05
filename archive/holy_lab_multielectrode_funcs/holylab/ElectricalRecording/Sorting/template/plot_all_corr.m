function plot_all_corr(spikeTimes, options)
%    if(length(spikeTimes)>5)
%       error('too many cells to evaluate');
%    end
   
   if(nargin<2)
      options=struct;
   end

   hfig=figure;
   nCells=length(spikeTimes);
   
   for row=1:nCells
      for col=1:row
        thisops = options;
        if isfield(options,'col')
          thisops.col = mean(options.col([row col],:),1);
        end
         hAxes=subplot(nCells, nCells, (row-1)*nCells+col);
         if(row==col)
           plot_corr(hAxes, spikeTimes{row},thisops);
         else
           plot_corr(hAxes, spikeTimes([row col]),thisops);
         end
         if(isfield(options, 'titles'))
            title(sprintf('(%s, %s)', options.titles{[row col]})); 
         end
      end % for, each col
   end % for, each row
   
   set(hfig, 'name', 'auto/cross correlation');
   