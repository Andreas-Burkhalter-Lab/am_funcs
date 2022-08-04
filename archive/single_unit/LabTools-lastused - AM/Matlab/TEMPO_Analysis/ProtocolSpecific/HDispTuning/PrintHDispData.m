%---------------------------------------------------------------------------------------------------------------------
%-- PrintHDispData.m : THis function prints out H Disp tuning curve specific data on the igure plot.
%--	BJP, 2/1/00
%---------------------------------------------------------------------------------------------------------------------

function   PrintHDispData(p_value, avg_resp, max_stats, min_stats, disp_groups, spont_level, cont_level, speeds, PATH, FILE, DDI, DTI, corr_coef, ASI);
axis([0 100 0 100]);
axis('off');
font_size = 8;
bump_size = 5.5;

for column = 1:size(speeds,1)           
    % type out stats onto screen
    xpos = -45 + (column)*35;   
    ypos = 20;
    line = sprintf('Speed: %g degrees/s', speeds(column));
    text(xpos,ypos,line,'FontSize',font_size+2);		ypos = ypos - bump_size;
    line = sprintf('MAX: Disp = %0.3g    Resp = %0.5g', max_stats(column).x, max_stats(column).y);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('MIN: Disp = %0.3g    Resp = %0.5g', min_stats(column).x, min_stats(column).y);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Avg. Resp. = %0.5g    Spont = %0.5g', avg_resp(column), spont_level);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Controls:  L = %0.3g   R = %0.3g   U = %0.3g', cont_level(column).left, cont_level(column).right, cont_level(column).uncorr);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Mod. Index (0 ->1) =  %0.5g', DTI(column));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('DDI =  %0.5g', DDI(column));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('ASI =  %0.5g', ASI(column));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('ANOVA: P = %0.3g', p_value(column));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size; 
    
    %print a summary line to the console, for batch use purposes
    buff = sprintf('%s %0.3g %0.4g %0.5g %0.4g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.3g %7.5f %0.4g\r\n', FILE, speeds(column), max_stats(column).x, ...
        max_stats(column).y, min_stats(column).x, min_stats(column).y, avg_resp(column), spont_level, ...
        cont_level(column).left, cont_level(column).right, cont_level(column).uncorr, DTI(column), DDI(column), corr_coef, p_value(column), ASI(column));
    disp(buff);
    
end

line = sprintf('Corr Coef: %0.3g', corr_coef);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size; 

return;