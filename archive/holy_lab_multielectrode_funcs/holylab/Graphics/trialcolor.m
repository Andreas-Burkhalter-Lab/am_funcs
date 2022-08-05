function [trialmap] = trialcolor(trialnumber, basemap)
%trialcolor produces an indexible cell array color map for use with multi-trial line plots where
%color is an indicator of trial order.  This function is just an algorithm for
%pulling out a subset of colors from a built-in color map at attractive intervals
%(determined by the number of trials in your data set). !!DON'T FORGET TO
%INDEX THE RETURNED COLORMAP AS A CELL ARRAY!!  See example below.
% Syntax:
%   [trialmap] = trialcolor(trialnumber, basemap)
% where
%   trialmap is the returned indexible colormap with a number of elements
%             matching the number of trials in your data set.
%   trialnumber is the number of trials (and therefore overlayed lines in your
%             multi-plot) from your experiment
%   basemap is an optional input that allows you to select the starting
%             color map from which the 'trialmap' will be calculated (default is set
%             to 'jet').
%usage example:  
%               [trialmap] = trialcolor(trialnumber, jet)
%               figure
%               hold on
%               for i = 1:trialnumber
%                   plot(xdata(i),ydata(i),'color',trialmap{i});
%               end
%               hold off
%
%Copywrite 2006 Terrence Holekamp
depth = 256;
switch lower(depth)
    case (trialnumber < 40)
        depth = 64;
    case trialnumber >= 40 && trialnumber < 73
        depth = 128;
end

if nargin < 2
    tempmap = jet(depth);
else
    tempmap = basemap(depth);
end

t = trialnumber;

if t > 128
    trialmap = tempmap;
else

    spread = round(length(tempmap)^(t/(t+3/t)));

    if t > 128
        n = 1;
    else
        n = round(spread/t);
    end

    trialmap = tempmap(1:n:spread,:);
end

rep = size(trialmap,1);
xslice = ones(1,rep);
trialmap =  mat2cell(trialmap, xslice, 3);




