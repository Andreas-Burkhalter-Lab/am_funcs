function scatcat(x,cat,factor,color,withmean,width)

% This script is going to make nice scatter plots with spaced points
% for discrete categories.

% x: vector with values to plot
% cat: category values for x-axis
% factor: spacing factor for the dots. The bigger the factor, the smaller
% the spacing (>10 approaches a straight line up and down).
% color: color of the dots

% Remove NaNs
k = find(isnan(x));
x(k) = [];
cat(k) = [];
l = length(x);


% This makes a l-sized vector c=[1 1 1 1 1 .... 1]

r = -1+2*rand(l,1); % random numbers to add to c
r = r./factor;
g = unique(cat);


xplace = cat+r;
m = length(color);

for i = 1:m
    scatter(xplace(cat==g(i)),x(cat==g(i)),'d','filled',color(i));
    hold on;
end

if nargin > 4
if withmean == 1
    [mn,s] = grpstats(x,cat,{'mean','sem'});
    if nargin > 5
    for i = 1:m
        semgraph(g(i),mn(i),s(i),'black',width);
    end
    else
    for i = 1:m
        semgraph(g(i),mn(i),s(i),'black',0.05);
    end  
    end
end
end

end


