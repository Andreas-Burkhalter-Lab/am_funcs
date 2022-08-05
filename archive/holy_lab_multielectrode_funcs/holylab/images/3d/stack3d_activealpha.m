function alpha = stack3d_activealpha(im,dfof,alphathresh,alphabg)
% stack3d_activealpha: pixelwise adjustment of alpha scaling based on dfof
%
% stack3d_activealpha allows you to set the relative opacity of active regions
%               of a 3-d image so that they can be visualized in depth
%               in the presence of non-active, but otherwise bright,
%               regions
% Examples of use:
%    alpha_func = @(im,dfof) stack3d_activealpha(im,dfof,alphathresh)
% or
%    alpha_func = @(im,dfof) stack3d_rgbalpha(im,dfof,alphathresh, alphabg)
%
% Where
%    alphathresh is a scalar between 0 and 1 indicating the 
%              threshold you would like to apply alpha changes 
%              (i.e., only pixels that have dfof > alphathresh
%              will have increased opacity over background 
%                  default value is 0.25
%    alphabg   is a scalar between 0 and 1 setting the alpha level of all
%              nonzero pixels which should be made "background" (i.e. transparent)
%                  default value is 0.25
%
%    NOTE: items must be supplied in this order, if you wish to fill a
%              later item (i.e. alphathresh or alphabg) but use the default
%              for others, supply the intermediates as []
%          For example, to use default alphathresh but alter alphabg to be 0.2,
%              you would type:
%                   alpha_func = @(im,dfof) stack3d_activealpha(dfof,[],0.2)
%
% The alpha_func would be supplied as an option to stack3d.
%
% See also: STACK3D.

% Copyright 2010 by Julian P. Meeks (Timothy Holy Laboratory)

if (nargin < 4)
    alphathresh = 0.25;
    alphabg = 0.25;
elseif (nargin < 3)
    alphabg = 0.25;
end

if isempty(dfof)
    dfof = zeros(size(im));
end

if isempty(alphathresh)
    alphathresh = 0.25;
end
if isempty(alphabg)
    alphabg = 0.25;
end
  
imkeep = im>nanmedian(nanmedian(nanmedian(im(im>0),3),2),1);

alpha = single(zeros(size(dfof)));
alpha(imkeep) = alphabg;
alpha(abs(dfof)>alphathresh) = alphabg+diff([alphabg 1])*abs(dfof(abs(dfof)>alphathresh));
alpha(~imkeep) = 0;
alpha(alpha < 0) = 0;
alpha(alpha > 1) = 1;
%alpha = alpha.*im;
end
