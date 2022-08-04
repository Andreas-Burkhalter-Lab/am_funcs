%---------------------------------------------------------------------------------------------------------------------
%-- PrintGratingContrastData.m : This function prints out Grating tuning curve specific data on the figure plot.
%--	CMA 06/05/06
%---------------------------------------------------------------------------------------------------------------------

function   PrintGratingContrastData2(p_value, Rmax, C50, n, DC_offset, max_rate, stats1, stats2, chi2, chiP, CDI, CMI);
axis([0 100 0 100]);
axis('off');
font_size = 9;
bump_size = 7;

xpos = -10;
ypos = 20;

line = sprintf('Fitted Tuning Curve:');
text(xpos,ypos,line,'FontSize',font_size+2,'Color','r');		ypos = ypos - bump_size;      
line = sprintf(' Rmax q(1) = %0.3g', Rmax);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;    
line = sprintf(' C50 q(2) = %0.3g', C50);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' n q(3) = %0.3g', n);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' DCoffset q(4) = %0.3g', DC_offset);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
line = sprintf(' Max = %0.3g', max_rate);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;         
% line = sprintf(' HWHM = %0.3g', sf_width_halfmax);
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       

xpos = 50;   
ypos = 10;
line = sprintf(' AnovP = %0.3g', p_value);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
line = sprintf(' CDI = %0.5g', CDI);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size; 
line = sprintf(' CMI = %0.5g', CMI);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' chi2 = %0.5g', chi2);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
line = sprintf(' chiP = %0.8g', chiP);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
line = sprintf(' RsqMn=%6.3f', stats1(:,1));
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' P=%8.6f', stats1(:,3));
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;

return;




% %---------------------------------------------------------------------------------------------------------------------
% %-- PrintGratingData.m : THis function prints out Grating tuning curve specific data on the figure plot.
% %--	GCD, 1/16/06
% %---------------------------------------------------------------------------------------------------------------------
% 
% function   PrintGratingData(p_value, avg_resp, max_stats, min_stats, disp_groups, spont_level, stim_types, PATH, FILE, DDI, DTI);
% axis([0 100 0 100]);
% axis('off');
% font_size = 8;
% bump_size = 5.5;
% 
% for column = 1:size(stim_types,1)           
%     % type out stats onto screen
%     xpos = -55 + (column)*55;   
%     ypos = 20;
%     line = sprintf('Gr. Type: %g (0=Sine,1=Sq)', stim_types(column));
%     text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
%     line = sprintf('MAX: Location = %0.3g    Resp = %0.5g', max_stats(column).x, max_stats(column).y);
%     text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
%     line = sprintf('MIN: Location = %0.3g    Resp = %0.5g', min_stats(column).x, min_stats(column).y);
%     text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
%     line = sprintf('Avg. Resp. = %0.5g    Spont = %0.5g', avg_resp(column), spont_level);
%     text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
%     line = sprintf('Mod. Index (0 ->1) =  %0.5g', DTI(column));
%     text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
%     line = sprintf('DDI =  %0.5g', DDI(column));
%     text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
%     line = sprintf('ANOVA: P = %0.3g', p_value(column));
%     text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size; 
%     
% end
% 
% return;