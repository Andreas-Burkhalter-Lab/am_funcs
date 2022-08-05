function msplotcompound(msdata,compound,options)
% MSPLOTCOMPOUND: show distribution of compound across fractions
% When you click on a point representing a fraction, it displays the mass
% spectrum for that fraction.  You can right-click on the plot to bring
% up a sliderwindow to zoom on the spectrum.
%
% Syntax:
%   msplotcompound(msdata,compound)
%   msplotcompound(msdata,compound,options)
% where
%   msdata is a structure of the type loaded my MSLOAD or MSCHOOSETRIALS;
%   compound is a struct of the type MSMATCHPEAK;
%   options is a structure with the following fields:
%     slope, offset: mass spec magnitudes are transformed according to
%           magout = slope*magin + offset
%       This is useful for aligned physiological data to ms data.
%
% See also: MSMATCHPEAK, MSLOAD, MSCHOOSETRIALS.

% Copyright 2005 by Timothy E. Holy

  slope = 1;
  offset = 0;
  if (nargin > 2)
    if isfield(options,'slope')
      slope = options.slope;
    end
    if isfield(options,'offset')
      offset = options.offset;
    end
  end

  % Plot the main window
  if isfield(compound,'stdmag')
    [xtmp,ytmp,errtmp] = addbreaks(compound.fraction,compound.meanmag,...
      compound.stdmag);
    hline = errorbar(xtmp,slope*ytmp+offset,slope*errtmp,'r-');
  else
    [xtmp,ytmp] = addbreaks(compound.fraction,compound.meanmag);
    hline = line(xtmp,slope*ytmp+offset,'Color','r');
  end
  set(hline,'HitTest','off');
  hold on
  % Set up the callbacks to display spectra when clicking on a data point
  nfraction = length(compound.fraction);
  for i = 1:nfraction
    hpoint = plot(msdata(i).fraction,slope*compound.meanmag(i)+offset,'ro');
    tmpdata = msdata(i);
    tmpdata.compoundlabel = compound.label;
    tmpdata.compoundmag = compound.meanmag(i);
    tmpdata.fraction = msdata(i).fraction;
    set(hpoint,'UserData',tmpdata,'MarkerFaceColor','r',...
               'MarkerSize',12,...
               'ButtonDownFcn',@mspc_plotspectrum);
  end
  hold off
  title(num2str(compound.label));
  xlabel('Fraction')
  ylabel('Ion current')
  set(gca,'XMinorGrid','on')
  
  
function mspc_plotspectrum(hObject,eventdata)
  msdata = get(hObject,'UserData');
  figure
  if (isfield(msdata,'isprof') && msdata.isprof)
    plot(msdata.rawmz(:,2),msdata.rawmz(:,1),'k');
  else
    stem(msdata.rawmz(:,2),msdata.rawmz(:,1),'k.');
  end
  hold on
  stem(msdata.compoundlabel,msdata.compoundmag,'r');
  hold off
  title(num2str(msdata.fraction));
  xlabel('m/z')
  ylabel('Ion current')
  slidwincmenu(gca);