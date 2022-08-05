% this will draw sems at the specific category x given plus and minus
% values and mn

function semgraph(x,mn,sem,color,width)

% So the user provides the following:
%           x: the user provided category. 
%           mn: the user provided mean
%           sem: the user provided sem
%           color: the color

% draw the mean:

line([(x-(2*width)) (x+(2*width))], [mn mn], 'Color', color, 'LineWidth',2.5);

% now draw the vertical line between plus and minus

hold on;

line([x x], [(mn - sem) (mn + sem)],'Color',color,'LineWidth',1.5);

% now draw the horizontal lines at the mn+sem and mn-sem

line([(x-width) (x+width)], [(mn+sem) (mn+sem)],'Color',color,'LineWidth',1.5);
line([(x-width) (x+width)], [(mn-sem) (mn-sem)],'Color',color,'LineWidth',1.5);




end
