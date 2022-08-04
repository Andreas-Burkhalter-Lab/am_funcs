%       This just adds the file name to the top of the plot.  RLS 10/05/07

function   PrintFileName(FILE);

axis([0 100 0 100]);
axis('off');
xpos = -10;
ypos = 110;
font_size = 8;
bump_size = 7;

temp = strcat(FILE);
temp(temp == '\') = '/';
% this prevents a stupid error from appearing on the screen
line = sprintf('File: %s', temp);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
return;