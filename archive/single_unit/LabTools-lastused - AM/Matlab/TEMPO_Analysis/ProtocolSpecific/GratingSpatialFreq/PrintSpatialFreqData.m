%The program changed to PrintSpatialFreqData.m to avoid issues with the
%protocol calling another function by this name.  Modified to print
%relevant data for Spatial Frequency protocol 04/28/06 CMA
%
%---------------------------------------------------------------------------------------------------------------------
%-- PrintGratingData.m : THis function prints out Grating tuning curve specific data on the figure plot.
%--	GCD, 1/16/06
%---------------------------------------------------------------------------------------------------------------------

function   PrintSpatialFreqData(p_value, base_rate, amplitude, peak, st_dev, log_offset, max_rate, stats1, stats2, chi2, chiP, SFDI, sf_width_halfmax);
axis([0 100 0 100]);
axis('off');
font_size = 9;
bump_size = 7;

xpos = -10;
ypos = 20;

line = sprintf('Fitted Tuning Curve:');
text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
line = sprintf(' BRate q(1) = %0.3g', base_rate);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;    
line = sprintf(' Ampl q(2) = %0.3g', amplitude);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf(' Peak q(3) = %0.3g', peak);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
line = sprintf(' SD q(4) = %0.3g', st_dev);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;         
line = sprintf(' Off q(5) = %0.3g', log_offset);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
line = sprintf(' HWHM = %0.3g', sf_width_halfmax);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       

xpos = 50;   
ypos = 10;
line = sprintf(' AnovP = %0.3g', p_value);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
line = sprintf(' SFDI = %0.5g', SFDI);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  
line = sprintf(' chi2 = %0.5g', chi2);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  
line = sprintf(' chiP = %0.8g', chiP);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  
line = sprintf(' RsqMn=%6.3f', stats1(:,1));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf(' P=%8.6f', stats1(:,3));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

return;