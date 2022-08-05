function ipo = imphysdf(ip,offset,weights, options)
% IMPHYSDF: calculate deltaf/f movies
%
% This function takes an input series of images "ip" and creates an
% output series "ipo" with the same number of images.  Each output image
% contains the normalized change in fluorescence.
%
% Syntax:
%   ipo = imphysdf(ip,offset,weight)
% where
%   ip is the input series of images (IMPHYS format)
%   offset is a cell array with two components: offset{1} specifies the
%     relative offset of images which contribute to the "foreground",
%     while offset{2} specifies the relative offset of images which
%     contribute to the "background"
%   weight is an optional cell array with two components: weight{1} is a
%     weights vector for offset{1}, and weight{2} for offset{2}.  By
%     default, the foreground images are averaged together, and the
%     background images are as well.  See IMPHYSTFILTER for a detailed
%     explanation of weighting.
% The output images are (foreground - background)./background.
%
% Example:
% Suppose you want to create a movie for which the foreground is defined
% by
%    fg(i) = (ip(i-1) + ip(i) + ip(i+1))/3
% and the background is defined by
%    bg(i) = (ip(i-9) + ip(i-10))/2.
%
% Then you'd call this function in the following way:
%    ipo = imphysdf(ip,{[-1 0 1],[-9 -10]})
%
% Frames near the beginning or end of the "movie" may not have
% the required predecessors or successors.  In that case, the
% closest frame is taken.  That is, a request for any frame < 1 will
% substitute frame 1, and any request for frame > "end" will use "end".
%
% VIMAGE note:
% In some cases you want to be able to identify an "ancestor" image. This
% is not a problem for functions that take only a single vimage argument,
% but it is obviously ambiguous for this function. The convention is that
% the _first_ vimage argument is the ancesctor.  Therefore, you might
% want to arrange your "offset" vector so that 0 appears first (and
% correspondingly rearrange "weight"), so that there is no
% temporal offset between the ancestor and the resultant
% image. Therefore, the above call might be better re-arranged as
%    ipo = imphysdf(ip,{[0 -1 1],[-9 -10]})
%
% See also: IMPHYSTFILTER, VIMAGE.
  
% Copyright 2005 by Timothy E. Holy

% An issue: should we use indexing, or stacknum? This currently uses
% indexing.
  
  nframes = length(ip);
  ipo = ip;    % Copy other fields
  if (nargin < 3 || (nargin>=3 && isempty(weights) ) )
    weights{1} = ones(1,length(offset{1}))/length(offset{1});
    weights{2} = ones(1,length(offset{2}))/length(offset{2});
  end
  if( nargin<4 )
     options=[]; 
  end
  
  if(~isfield(options, 'norm'))
     options.norm=1; 
  end
  
  for i = 1:nframes
    indexfg = ipdf_boundaries(i + offset{1},nframes);
    indexbg = ipdf_boundaries(i + offset{2},nframes);
    framesfg = {ip(indexfg).image};
    framesbg = {ip(indexbg).image};
    fg = imwsum(weights{1},framesfg{:});
    bg = imwsum(weights{2},framesbg{:});
    if(options.norm)
       ipo(i).image = (fg-bg)./bg;
       ipo(i).computation='\DeltaF/F';
    else
       ipo(i).image = fg-bg;
       ipo(i).computation='\DeltaF';
    end
  end
  
function index = ipdf_boundaries(index,nframes)  
  index = max([index; ones(1,length(index))]);  % None < 1
  index = min([index; nframes*ones(1,length(index))]); % None > nframes
