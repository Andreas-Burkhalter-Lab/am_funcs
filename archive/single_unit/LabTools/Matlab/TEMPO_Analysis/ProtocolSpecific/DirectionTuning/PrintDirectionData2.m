%---------------------------------------------------------------------------------------------------------------------
%-- PrintDirectionData.m : THis function prints out Direction tuning curve specific data on the igure plot.
%--	BJP, 2/2/00
%---------------------------------------------------------------------------------------------------------------------

function   PrintDirectionData2(p_value, base_rate, spont_level, amplitude, pref_dir, max_rate, width, DSI, stats1, stats2, DirDI, chi2, chiP)
axis([0 100 0 100]);
axis('off');
font_size = 9;
bump_size = 7;

xpos = -10;   
ypos = 20;

line = sprintf('Fitted Tuning Curve:');
text(xpos,ypos,line,'FontSize',font_size+2,'Color','r');		ypos = ypos - bump_size;      
line = sprintf('Preferred Direction = %1.5g   Resp = %0.5g', pref_dir, max_rate);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;    
line = sprintf('Base Rate = %0.5g   Spont = %0.5g', base_rate, spont_level);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf('Amplitude = %0.3g', amplitude);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
line = sprintf('Width (FWHM) = %0.3g', width);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;         
line = sprintf('DSI = %0.3g', DSI);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       

xpos = 60;   
ypos = 10;
font_size = 8;
line = sprintf('ANOVA: P = %0.3g', p_value);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
line = sprintf('Dir Discrim Ind = %0.5g', DirDI);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
line = sprintf('chi2 = %0.5g   chiP = %0.8g', chi2, chiP);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1(1), stats1(3));
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2(1), stats2(3));
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;

return;