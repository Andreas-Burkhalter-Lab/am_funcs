%---------------------------------------------------------------------------------------------------------------------
%-- PrintGratingData.m : THis function prints out Grating Temporal Frequency tuning curve specific data on the igure plot.
%--	CMA 04/10/06
%---------------------------------------------------------------------------------------------------------------------

function PrintSpeedData2(p_value, base_rate, amplitude, peak, st_dev, log_offset, max_rate, stats1, stats2, chi2, chiP, SDI, HWHMR, HWHML);
axis([0 100 0 100]);
axis('off');
font_size = 9;
bump_size = 7;

xpos = -10;
ypos = 20;

line = sprintf('Fitted Tuning Curve:');
text(xpos,ypos,line,'FontSize',font_size+2,'Color','r');		ypos = ypos - bump_size;      
line = sprintf(' BRate q(1) = %0.3g', base_rate);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;    
line = sprintf(' Ampl q(2) = %0.3g', amplitude);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' Peak q(3) = %0.3g', peak);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
line = sprintf(' SD q(4) = %0.3g', st_dev);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;         
line = sprintf(' Off q(5) = %0.3g', log_offset);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
line = sprintf(' HWHML = %5.3f', HWHML);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
line = sprintf(' HWHMR = %5.3f', HWHMR);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       


xpos = 50;   
ypos = 10;
line = sprintf(' AnovP = %0.4f', p_value);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
line = sprintf(' TFDI = %0.5g', SDI);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
line = sprintf(' chi2 = %0.5g', chi2);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
line = sprintf(' chiP = %0.8g', chiP);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
line = sprintf(' RsqMn=%6.3f', stats1(1));
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' P=%8.6f', stats1(3));
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;

return;




% %---------------------------------------------------------------------------------------------------------------------
% %-- PrintSpeedData.m : THis function prints out Speed tuning curve specific data on the igure plot.
% %--	BJP, 2/2/00
% %---------------------------------------------------------------------------------------------------------------------
% 
% function   PrintSpeedData(speeds, spk_rates, max_stats, min_stats, speed_groups, spont_level, fit_params, p_value, avg_resp, stats1, stats2, SDI, chi2, chiP);
% 
% axis([0 100 0 100]);
% axis('off');
% font_size = 9;
% bump_size = 7;
% 
% % type out stats
% xpos = -10;   
% ypos = 20;
% line = sprintf('Fitted Tuning Curve:');
% text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
% line = sprintf('MAX: Speed = %1.3g   Resp = %0.5g', max_stats.x, max_stats.y);
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
% line = sprintf('MIN: Speed = %1.3g   Resp = %0.5g', min_stats.x, min_stats.y);
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
% line = sprintf('Avg. Resp. = %0.5g   Spont = %0.5g', avg_resp, spont_level);
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
% line = sprintf('Tuning ANOVA: P = %0.3g', p_value);
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
% line = sprintf('Spd Discrim Ind = %0.5g', SDI);
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
% line = sprintf('chi2 = %0.5g   chiP = %0.8g', chi2, chiP);
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
% 
% font_size = 8;
% bump_size = 6.5;
% xpos = 55;   
% ypos = 30;
% line = sprintf('Fitting Parameters:');
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
% line = sprintf('Base Rate q(1) = %0.3g', fit_params(1));
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;        
% line = sprintf('Amplitude q(2) = %0.3g', fit_params(2));
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
% line = sprintf('Exp Rise q(3) = %0.3g', fit_params(3));
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
% line = sprintf('Offset q(4) = %0.3g', fit_params(4));
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
% line = sprintf('Exp Decay q(5) = %0.3g', fit_params(5));
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;          
% % line = sprintf('Fitting Error: %0.3g', fit_error);
% % text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
% line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1(1), stats1(3));
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
% line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2(1), stats2(3));
% text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
% 
% return;