%---------------------------------------------------------------------------------------------------------------------
%-- PrintSpeedData.m : THis function prints out Speed tuning curve specific data on the igure plot.
%--	BJP, 2/2/00
%---------------------------------------------------------------------------------------------------------------------

function   PrintSpeedDataLogGauss(speeds, spk_rates, max_stats, min_stats, speed_groups, spont_level, fit_params, fit_error, p_value, avg_resp, stats1, stats2, SDI, chi2, chiP);

axis([0 100 0 100]);
axis('off');
font_size = 9;
bump_size = 7;

% type out stats
xpos = -10;   
ypos = 20;
line = sprintf('Fitted Tuning Curve:');
text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;      
line = sprintf('MAX: Speed = %1.3g   Resp = %0.5g', max_stats.x, max_stats.y);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('MIN: Speed = %1.3g   Resp = %0.5g', min_stats.x, min_stats.y);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Avg. Resp. = %0.5g   Spont = %0.5g', avg_resp, spont_level);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Tuning ANOVA: P = %0.3g', p_value);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  
line = sprintf('Spd Discrim Ind = %0.5g', SDI);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  
line = sprintf('chi2 = %0.5g   chiP = %0.8g', chi2, chiP);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  

font_size = 8;
bump_size = 6.5;
xpos = 55;   
ypos = 30;
line = sprintf('Fitting Parameters:');
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
line = sprintf('Base Rate q(1) = %0.3g', fit_params(1));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;        
line = sprintf('Amplitude q(2) = %0.3g', fit_params(2));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;  
line = sprintf('Peak Speed q(3) = %0.3g', fit_params(3));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
line = sprintf('Sigma q(4) = %0.3g', fit_params(4));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
line = sprintf('Log Offset q(5) = %0.3g', fit_params(5));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;          
line = sprintf('Fitting Error: %0.3g', fit_error);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;       
line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1(1), stats1(3));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2(1), stats2(3));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

return;