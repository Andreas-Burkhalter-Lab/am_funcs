%---------------------------------------------------------------------------------------------------------------------
%-- PrintGratingOrientationData.m : THis function prints out Grating tuning curve specific data on the figure plot.
%--	GCD, 1/16/06
%---------------------------------------------------------------------------------------------------------------------

function   PrintGratingOrientationData2(p_value, avg_resp, max_stats, min_stats, disp_groups, spont_level, stim_types, PATH, FILE, DDI, HWHM);
axis([0 100 0 100]);
axis('off');
font_size = 8;
bump_size = 7;

for column = 1:size(stim_types,1)           
    % type out stats onto screen
    xpos = -65 + (column)*55;   
    ypos = 20;
    line = sprintf('Gr. Type: %g (0=Sine,1=Sq)', stim_types(column));
    text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
    line = sprintf('MAX: Location = %0.3g    Resp = %0.5g', max_stats{column}.x, max_stats{column}.y);
    text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
    line = sprintf('MIN: Location = %0.3g    Resp = %0.5g', min_stats{column}.x, min_stats{column}.y);
    text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
    line = sprintf('Avg. Resp. = %0.5g    Spont = %0.5g', avg_resp(column), spont_level);
    text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
    %line = sprintf('Mod. Index (0 ->1) =  %0.5g', DTI(column));
    %text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('DDI =  %0.5g', DDI(column));
    text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
    line = sprintf('ANOVA: P = %0.3g', p_value(column));
    text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size; 
    line = sprintf('HWHM = %0.5g', HWHM);
    text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size; 
    
end

return;