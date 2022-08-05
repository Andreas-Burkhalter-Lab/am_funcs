function drgui(ephysstim,fieldname,trange,valvelabels)
% DRGUI: a graphical browser for electrophysiology results
% This creates an image in which individual pixels correspond to
% cell/stimulus (or electrode/stimulus) pairs; clicking on a pixel brings
% up an ephysgui browser so you can look at PSTHs, rasters, etc.
%
% Syntax:
%   drgui(ephysstim,fieldname,trange)
%   drgui(ephysstim,fieldname,trange,valvelabels)
% where
%   ephysstim is your array of ephys structures, tagged by stimulus
%   fieldname is 'sniptimes' or 'celltimes', depending on whether you want
%     to look at multiunit data ('sniptimes') or single units ('celltimes')
%   trange is a 3-vector such as [-5 0 20], which specifies the deltar time
%     intervals (see DELTARATE)
%   valvelabels is an optional cell array of stimulus names; by default,
%     all stimuli will be used.  One reason to supply this argument is to
%     control the order of presentation of stimuli on the y-axis.
%
% See also: DELTARATEGUI.

% Copyright 2007 by Timothy E. Holy
% July 2011 changed firing rate metric to drmax_across trials (HAA)

  if (nargin < 4)
    valvelabels = unique({ephysstim.tag});
  end
  
  % Create the proper column data for mdexplore_im_ephys
  coldata = struct;
  switch(fieldname)
   case 'celltimes'
     % Set it up so that clicking on a point shows you the cell, but have
     % it default to the correct channel if the user wants to see the
     % sniptimes
%      for i = 1:length(ephysstim(1).cellnums)
%        coldata(i).field = {'cellnumber','channelnumber'};
%        coldata(i).value = {ephysstim(1).cellnums(i),ephysstim(1).cellchandef{i}(1)};
%      end
     % Now: ephysgui will figure out which channel to use if you don't
     % supply it. That's better, because now the channel will follow
     % changes in the cellnumber.
     for i = 1:length(ephysstim(1).cellnums)
       coldata(i).field = 'cellnumber';
       coldata(i).value = ephysstim(1).cellnums(i);
     end

   case 'sniptimes'
     % Set it up so that clicking on a point shows you the channel, but if
     % cells are available, then have it default to the first cell on that
     % channel (or the first cell, if no cell is defined on that channel)
     for i = 1:length(ephysstim(1).channels)
       this_channel = ephysstim(1).channels(i);
       if isfield(ephysstim(1),'celltimes')
         cellIndex = find(round(ephysstim(1).cellnums) == this_channel);
         if ~isempty(cellIndex)
           cellnum = ephysstim(1).cellnums(cellIndex(1));
         else
           cellnum = ephysstim(1).cellnums(1);
         end
         coldata(i).field = {'channelnumber','cellnumber'};
         coldata(i).value = {this_channel,cellnum};
       else
         coldata(i).field = 'channelnumber';
         coldata(i).value = this_channel;
       end
     end
   otherwise
    error(['Fieldname ' fieldname ' not recognized']);
  end
  % Create the row data (info about stimuli)
  rowdata = struct('field','tags','value',valvelabels);

  % Set up the plotting function
  pp = struct('fieldtoplot',fieldname,'type','PSTH w/ sem','showtags', ...
	      true,'raster_rmax_times',[0 1]);
  params.plotfunc = @(col,row) ...
      mdexplore_im_ephys(col,row,ephysstim,pp,coldata,rowdata);

  % Now we have to calculate the deltars to show in the image
  % [dr,drerr,tmp,npb,n_trials] = deltarate(ephysstim,valvelabels,trange,fieldname);
  [dr,drerr,tmp,npb,n_trials,tmax] = drmax_acrosstrials(ephysstim,valvelabels,trange,fieldname);
  % Create a normalized response, as dr/max(abs(dr)+drerr) for each cell
%   drdenom = max(max(abs(dr),abs(dr) + drerr.*sqrt(n_trials)),[],1);
  mserr = nanmean(drerr.^2,1);
  drerr_regularized = sqrt(bsxfun(@plus,drerr.^2,mserr));
  drdenom = max(abs(dr)+drerr_regularized,[],1);
  % Ones for which drerr is undefined just get normalized by their max rate
  disnan = isnan(drdenom);
  drdenom(disnan) = max(abs(dr(:,disnan)),[],1); 
  % Ones that never spiked should be sent to 0
  drdenom(drdenom == 0) = 1;
  dr(isnan(dr)) = 0;
  drnorm = dr ./ repmat(drdenom,size(dr,1),1);
  
  % Do the plot
  params.ylabel = valvelabels;
  figure;
  hax = SplitVert([0.9 0.92],[1 0 1]);
  params.hax = hax(2);
  plot(hax(1),max(dr./drerr_regularized,[],1));
  axis(hax(1),'tight');
  mdexplore_im(drnorm,params);
  set(hax(1),'XLim',get(params.hax,'XLim'));
  
  % Set up the colormap
  colormap(params.hax,colormap_pm(drnorm));
  colorbar('peer',params.hax)
  posim = get(params.hax,'Position');
  posplot = get(hax(1),'Position');
  posplot(3) = posim(3);
  set(hax(1),'Position',posplot,'XTick',[]);
  
  
  