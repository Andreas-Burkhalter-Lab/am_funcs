function ret = imphysplot(ip,pp,hax)
% IMPHYSPLOT: show data plots for imaging analysis
% Syntax:
%   ret = imphysplot(ip,pp)
%   ret = imphysplot(ip,pp,hax)
% where
%   ip is the input imphys array;
%   pp is the input IMPHYSPLOTPARAMS structure;
%   hax (default current axis) is the handle of the axis to use for
%     plotting;
% and
%   ret is returned data, the identity depending on the plot mode (see
%     IMPHYSPLOTPARAMS). 
%
% See IMPHYSPLOTPARAMS for details about the types of plots and how you
% specify their parameters.
%
% See also: IMPHYSPLOTPARAMS.
  
% Copyright 2005 by Timothy E. Holy
  
  % Axis preliminaries
  if (nargin < 3)
    hax = gca;
  end
  % Check input axis handle
  if (~ishandle(hax) | ~strcmp(get(hax,'Type'),'axes'))
    error('Input axis is not valid');
  end
  % Delete any legends that are hanging around in this axis
  hleg = findobj(0,'Tag','legend');
  for i = 1:length(hleg)
    ud = get(hleg(i),'UserData');
    if (ud.PlotHandle == hax)
      delete(hleg(i));
    end
  end
  % This is a cla for axis hax
  delete(findobj(get(hax,'Children'),'flat','HandleVisibility','on'))
  
  % Convert mode to lower case
  pp.mode = lower(pp.mode);
  
  % Process different modes
  switch pp.mode
    case {'image','image w/ roi'},
     % Compute the background, if desired
     if (isfield(pp,'bgindx') & ~isempty(pp.bgindx))
       imbg = double(imphysfetch(ip(pp.bgindx(1))));
       for i = 2:length(pp.bgindx)
         imbg = imbg + double(imphysfetch(ip(pp.bgindx(i))));
       end
       imbg = imbg/length(pp.bgindx);   % Mean value
     else
       imbg = [];
     end
     % Compute the foreground (average over frames)
     imfg = double(imphysfetch(ip(pp.fgindx(1))));
     for i = 2:length(pp.fgindx)
       imfg = imfg + double(imphysfetch(ip(pp.fgindx(i))));
     end
     imfg = imfg/length(pp.fgindx);   % Mean value
     % Do the subtraction
     if ~isempty(imbg)
       im = imfg - imbg;
     else
       im = imfg;
     end
     % Do spatial filtering (do this before cropping to get edges right,
     % even if it's slower)
     if isfield(pp,'spatialfilter')
       im = imfilter(im,pp.spatialfilter);
     end
     % Set up a temporary ephys structure to handle dimensioning issues
     imip = imphyscopy(ip(pp.fgindx(1)),'dimension');  % copy size params 
     imip.image = im;
     % Crop the image (use imphyscrop for consistency)
     if isfield(pp,'cropstruct')
       imip = imphyscrop(imip,pp.cropstruct);
     end
     % Plot the image
     if isfield(pp,'clim')
       imagesc(imip.xrange,imip.yrange,im,clim,'parent',hax);
     else
       imagesc(imip.xrange,imip.yrange,im,'parent',hax);
     end
     % Set up the return value
     if (nargout > 0)
       ret = im;
     end
     
     % If plotting ROIs, get on with it
     if strcmp(pp.mode,'image w/ roi')
       if isfield(pp,'roioptions')
         plotroi(pp.roidata,pp.roioptions);
       else
         plotroi(pp.roidata);
       end
     end
     
     % Tweak appearance
     set(gca,'TickDir','out');
     
   case 'stimulus',
    % Use the ephys plotting routines
    stim = [ip.stimulus; ip.stacktime];
    cstimindx = unique([1, find(diff([ip.stimulus]))+1, size(stim,2)]);
    ep.stimulus = stim(:,cstimindx);
    if ~isfield(pp,'type')
      pp.type = 'bar';
    end
    epp = struct('fieldtoplot','stimulus','type',pp.type);
    ephysplot(ep,epp);
   
   case 'psth',
    error('Not implemented yet')
  end
  
