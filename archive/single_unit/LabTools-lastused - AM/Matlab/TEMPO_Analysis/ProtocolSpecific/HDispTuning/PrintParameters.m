function PrintParameters(data, pars, stats1, stats2, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, freq, p_value, avg_resp, max_stats, min_stats, disp_groups, spont_level, cont_level, speeds, DDI, DTI, corr_coef, ASI);

	ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

	axis([0 100 0 100]);
	axis('off');
   font_size = 9;
   bump_size = 5;
   
   xpos = 55;
   ypos = 20;
   
	for column = 1:size(speeds,1)           
   	% type out stats onto screen
		xpos = 40 + (column)*25;   
   	ypos = 20;
   line = sprintf('base rate: %f',pars{column}(1));
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
   line = sprintf('amplitude: %f',pars{column}(2));
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
   line = sprintf('center: %f',pars{column}(3));
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
   line = sprintf('size: %f',pars{column}(4));
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
   line = sprintf('frequency: %f',pars{column}(5));
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
   line = sprintf('phase: %f',pars{column}(6));
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
   line = sprintf('FT frequency: %f',freq(column));
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - 2*bump_size;

	line = sprintf('R-squared mean: %f', stats1{column}(1));
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
   line = sprintf('R-squared raw: %f', stats2{column}(1));
   text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
   
   fid = fopen('Z:\Data\Tempo\Batch Files\GaborFits.out','a');
   %prints base rate, amplitude, center, size, frequency, phase, FT frequency, R-squared mean, R-squared raw, pref. direction, pref. speed, RF X-Ctr, RF Y-Ctr, RF diameter, Spont, Controls (L,R,U), Mod.Index(0->1), DDI, ASI, ANOVA
   fprintf(fid,'%s %7.2f %7.2f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %5.2f %5.2f %5.2f %5.2f %5.2f %8.5f %8.3f %8.3f %8.3f %8.5f %8.5f %8.5f %8.5f\n',FILE, pars{column}(1), pars{column}(2), pars{column}(3), pars{column}(4), pars{column}(5), pars{column}(6), freq(column), stats1{column}(1), stats2{column}(1), data.one_time_params(PREFERRED_DIRECTION), speeds(column), data.one_time_params(RF_XCTR), data.one_time_params(RF_YCTR), data.one_time_params(RF_DIAMETER), spont_level, cont_level(column).left, cont_level(column).right, cont_level(column).uncorr, DTI(column), DDI(column), ASI(column), p_value(column));
   fclose(fid);
   
   end
    
   
return;