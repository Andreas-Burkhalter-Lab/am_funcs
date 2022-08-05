function inspect_tube_traversals(intensityfilename,options)
% inspect_tube_traversals: check movie correspondence to automated behavioral data
%
% Allow a user to see what was happening in the video by clicking at
% certain points in the timecourse.
%
% Syntax:
%   inspect_tube_traversals
%   inspect_tube_traversals(intensityfilename)
%   inspect_tube_traversals(intensityfilename,options)
% where
%   intensityfilename is a string containing the name of the _intensity.mat
%     file output by track_tube_traversals. If not supplied, or supplied as
%     empty, the user is prompted to pick the file from the file chooser.
%   options is a structure which may have the following fields:
%     raw (default false): if true, shows the raw intensity ratios, rather
%       than the parsed entrances and exits.
%
% See also: track_tube_traversals.

% Copyright 2010 by Timothy E. Holy

  if (nargin < 1 || isempty(intensityfilename))
    intensityfilename = uigetfile('*_intensity.mat');
  end
  if (nargin < 2)
    options = struct;
  end
  options = default(options,'raw',false);
  s = load(intensityfilename);
  mov = mmreader_ffmpeg(s.videofile,'calibrate',true);
  figure;
  if options.raw
    n_circles = size(s.circles,1)/2;
    rat = s.intensity(:,1:n_circles,:) ./ s.intensity(:,n_circles+1:end,:);
    hline = plot(s.trange(1):s.trange(2),rat(:,[1 2],1));
  else
    [loc,t] = analyze_tube_traversals(s.intensity);
    hline = stairs(t+s.trange(1),loc);
  end
  controls = mplay(mov);
  set(hline,'HitTest','off')
  set(gca,'YLim',[-0.2 2.2]);
  sliderwindow(gca)
  set(gca,'ButtonDownFcn',@(src,evt) itt_changeframe(src,evt,controls));
  
  function itt_changeframe(src,~,controls)
    cp = get(src,'CurrentPoint');
    s = struct(controls);
    s.fcns.jumptoframe(s.hfig,round(cp(1)));
%     controls.jumptoframe(cp(1));
    
