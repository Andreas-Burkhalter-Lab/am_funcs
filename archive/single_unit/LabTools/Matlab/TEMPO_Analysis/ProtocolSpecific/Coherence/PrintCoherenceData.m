%The program changed to PrintContrastData.m to avoid issues with the
%protocol calling another function by this name.
%
%---------------------------------------------------------------------------------------------------------------------
%-- PrintGratingData.m : THis function prints out Grating tuning curve specific data on the figure plot.
%--	GCD, 1/16/06
%-- Modified for coherence protocol on 11/27/06 by CMA
%---------------------------------------------------------------------------------------------------------------------

function   PrintCoherenceData(unique_direc, p, Rsq_means, P_means, Rsq_raw, P_raw)
%(p_value, base_rate, amplitude, peak, st_dev, log_offset, max_rate, stats1, stats2, chi2, chiP, SFDI, sf_width_halfmax);
axis([0 100 0 100]);
axis('off');
font_size = 9;
bump_size = 7;

xpos = -10;
ypos = 20;

line = sprintf('Fitted Tuning Curve:');
text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
line = sprintf(' Directions = %0.3g', unique_direc);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;    
line = sprintf(' Slope = %0.3g', p(:,1));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf(' Intercept = %0.3g', p(:,2));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       

xpos = 50;   
ypos = 10;
line = sprintf(' Rsqmeans = %0.3g', Rsq_means);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
line = sprintf(' Pmeans = %0.5g', P_means);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  
line = sprintf(' Rsqraw = %0.5g', Rsq_raw);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  
line = sprintf(' Praw = %0.8g', P_raw);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  


return;