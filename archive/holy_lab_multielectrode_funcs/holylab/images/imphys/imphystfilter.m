function ipo = imphystfilter(ip,offset,weight)
% IMPHYSTFILTER: temporally filter a series of images
%
% This function takes an input series of images "ip" and creates an
% output series "ipo" with the same number of images.  Each output image
% is a linear combination of input images.
%
% Syntax:
%   ipo = imphystfilter(ip,offset,weight)
% where
%   ip is the input series of images (IMPHYS format)
%   offset specifies the relative offset of images which contribute to
%     the output;
%   weight specifies how the images refered to by "offset" are to be
%     weighted.
%
% Example:
% Suppose you want to create a movie for which
%    ipo(i) = (ip(i-1) + ip(i) + ip(i+1))/3 - (ip(i-5) + ip(i-4))/2
% (that is, the output is the average of three adjacent frames, minus the
% average of two earlier frames).  Then you'd call this function with
%   offset =  [-5    -4    -1    0    1]
%   weight =  [-0.5  -0.5  1/3   1/3  1/3]
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
% temporal offset between the ancestor and the resultant image.
%
% See also: VIMAGE
  
% Copyright 2005 by Timothy E. Holy

% An issue: should we use indexing, or stacknum? This currently uses
% indexing.
  
  nframes = length(ip);
  ipo = ip;    % Copy other fields
  for i = 1:nframes
    index = i + offset;
    index = max([index; ones(1,length(offset))]);
    index = min([index; nframes*ones(1,length(offset))]);
    frames = {ip(index).image};
    ipo(i).image = imwsum(weight,frames{:});
  end
  