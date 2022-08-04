function PrintDataArray(variable_data, variable_names)
%------------------------------------------------------------------------
% this function takes two cell arrays.  one cell array consists of arrays 
% of the data that will be printed to the screen.  the other cell array
% consists of the corresponding labels to be printed with the data.
% this function will then take both cell arrays and print them into columns
% onto a specified figure.  2/8/00 - JDN
%------------------------------------------------------------------------

axis([0 100 0 100]);
axis('off');
xpos = -10;
startypos = 110;
font_size = 9;
bump_size = 10;
   
for i = 1:length(variable_names)
	xpos = xpos+15;
   ypos = startypos;
   text(xpos, ypos, char(variable_names(i)), 'FontSize', font_size);
   temp = variable_data{i};
   for j = 1:length(temp)
   	ypos = ypos - bump_size;
      line = sprintf('%3.3f', temp(j));
		text(xpos,ypos,line,'FontSize',font_size);
   end
end