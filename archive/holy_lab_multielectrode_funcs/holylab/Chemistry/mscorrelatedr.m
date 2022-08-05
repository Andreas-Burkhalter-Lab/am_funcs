function mscorrelatedr(msdata,compound,drstruct,options)
% MSCORRELATEDR: correlate neural activity with msdata
% Syntax:
%   mscorrelatedr(msdata,compound,drstruct)
% where
%   msdata is a structure array containing mass spectrum data of the type
%     loaded by MSLOAD or MSCHOOSETRIALS;
%   compound is a structure array containing information about individual
%     compounds across fractions, of the type returned by MSGETCOMPOUNDS;
%   drstruct is a structure containing delta-rate data, with the
%     following fields:
%       fraction: a vector of fraction numbers
%       channel: a vector of channel numbers/indices
%       dr: a nfractions-by-nchannels matrix containing deltar for each
%         fraction/channel pair;
%       drerr: a nfractions-by-nchannels matrix containing the standard
%         errors in dr for each fraction/channel pair.
%   options is a structure which may have the following fields:
%       userank: if true, correlation is performed on the basis of ranks,
%         rather than actual values (default: false)
%
% When there are multiple compounds and multiple channels, the
% correlations are presented as a matrix, where color indicates the
% degree of correlations.  When the input instead has either a single
% compound OR a single channel, then the correlation is presented as a
% graph (making it easier to see the quantitative aspects).
  
% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 4)
      options = struct;
  end
  if ~isfield(options,'userank')
      options.userank = 0;
  end
  
  % Find common fractions
  [cfrac,idr,icomp] = intersect(drstruct.fraction,compound(1).fraction);
  % Correlate activity of channels with compound concentrations
  drtmp = drstruct.dr(idr,:);
  if options.userank
    % Try rank correlation
    for i = 1:size(drtmp,2)
      drtmp(:,i) = midrank(drtmp(:,i));
    end
  end
  for i = 1:length(compound)
    tcompound = compound(i).meanmag(icomp);
    %if options.userank
    %    tcompound = midrank(tcompound);  % Don't do this, because MS data
    %                                     % do not saturate
    %end
    [slope(i,:),offset,r(i,:)] = ...
        linregress(repmat(tcompound(:),1,size(drtmp,2)),drtmp);
  end
  % Store data in figure
  hfig = figure;
  setappdata(hfig,'msdata',msdata);
  setappdata(hfig,'compound',compound);
  setappdata(hfig,'drstruct',drstruct);
  setappdata(hfig,'slope',slope);
  
  if (length(compound) == 1)
    % Compound-centric: pick a compound, look for channels whose activity
    % correlates with the compound's concentration
    hold on
    for i = 1:length(r)
      hpoint = plot(i,r(i),'bo');
      set(hpoint,'MarkerFaceColor','b',...
                 'MarkerSize',12,...
                 'UserData',[i 1],...
                 'ButtonDownFcn',@mscd_handler1d);
    end
    hold off
    xlabel('Channel index');
    ylabel('Correlation with compound');
    title(['Compound ' num2str(compound.label)])
  elseif (length(drstruct.channel) == 1)
    % Channel-centric: pick a channel, look for compounds whose
    % concentration correlates with the channel's response
    hold on
    for i = 1:length(compound)
      hpoint = plot(i,r(i),'bo');
      set(hpoint,'MarkerFaceColor','b',...
                 'MarkerSize',12,...
                 'UserData',[1 i],...
                 'ButtonDownFcn',@mscd_handler1d);
    end
    hold off
    xlabel('Compound index');
    ylabel('Correlation with \Deltar');
    title(['Channel ' num2str(drstruct.channel)])
  else
    % Color plot of both channels and compounds
    imagesc(r,[-1,1]);
    colorbar
    xlabel('Channel index');
    ylabel('Compound index');
    install_mouse_event_handler(gca, 'down', @mscd_handler2d);
  end

  
function mscd_plotfractions(hfigparent,pos)
  % Note pos(1) = channel index, pos(2) = compound index
  msdata = getappdata(hfigparent,'msdata');
  compound = getappdata(hfigparent,'compound');
  drdata = getappdata(hfigparent,'drstruct');
  slope = getappdata(hfigparent,'slope');
  figure
  [xtmp,ytmp,etmp] = addbreaks(drdata.fraction,drdata.dr(:,pos(1)),...
                               drdata.drerr(:,pos(1)));
  errorbar(xtmp,ytmp,etmp,'k');
  hold on
  %msplotcompound(msdata,compound(pos(2)),...
  %               struct('slope',abs(slope(pos(2),pos(1)))));
  msplotcompound(msdata,compound(pos(2)),...
                 struct('slope',slope(pos(2),pos(1))));
  title(['Channel ' num2str(drdata.channel(pos(1))) ...
         ', compound ' num2str(compound(pos(2)).label)]);
  ylabel('\Deltar')
  
function mscd_handler1d(hObject,eventdata)
  pos = get(hObject,'UserData');
  mscd_plotfractions(gcf,pos);
  
function ret = mscd_handler2d(src,param)
  pos = get(gca,'currentpoint');
  pos = round(pos(1,1:2));
  mscd_plotfractions(gcf,pos);
  ret = 1;
  
