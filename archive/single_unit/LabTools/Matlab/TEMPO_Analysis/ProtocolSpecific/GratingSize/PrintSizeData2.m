%---------------------------------------------------------------------------------------------------------------------
%-- PrintSizeData.m : This function prints out Grating tuning curve specific data on the figure plot.
%--	CMA 06/12/06
%---------------------------------------------------------------------------------------------------------------------

function   PrintSizeData2(p_value, Ke, a, Ki, b_plus_a, DC, stats_mean, stats_raw, OptSize, SI, SDI, Erf_chi2, Erf_chiP, Diff_Erf_chi2, Diff_Erf_chiP);
axis([0 100 0 100]);
axis('off');
font_size = 9;
bump_size = 7;

xpos = -10;
ypos = 25;

line = sprintf('Fitted Tuning Curve:');
text(xpos,ypos,line,'FontSize',font_size+2,'Color','r');		ypos = ypos - bump_size;      
line = sprintf(' Ke q(1) = %0.3g', Ke);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;    
line = sprintf(' a q(2) = %0.3g', a);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' Ki q(3) = %0.3g', Ki);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
line = sprintf(' b+a q(2)+q(4) = %0.3g', b_plus_a);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;         
line = sprintf(' DCOff q(5) = %0.3g', DC);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
line = sprintf(' OptSize = %0.3g', OptSize);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' SI=%6.3f', SI);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' AnovP = %0.3g', p_value);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;       
 

xpos = 50;   
ypos = 25;

line = sprintf(' SDI = %0.5g', SDI);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size; 
line = sprintf(' Erf chi2 = %0.5g', Erf_chi2);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
line = sprintf(' Erf chiP = %0.8g', Erf_chiP);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
line = sprintf(' DiffErf chi2 = %0.5g', Diff_Erf_chi2);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
line = sprintf(' DiffErf chiP = %0.8g', Diff_Erf_chiP);
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;  
line = sprintf(' RsqMn=%6.3f', stats_mean(:,1));
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' PMn=%8.6f', stats_mean(:,3));
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' RsqRaw=%6.3f', stats_raw(:,1));
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;
line = sprintf(' PRaw=%8.6f', stats_raw(:,3));
text(xpos,ypos,line,'FontSize',font_size,'Color','r');		ypos = ypos - bump_size;


return;