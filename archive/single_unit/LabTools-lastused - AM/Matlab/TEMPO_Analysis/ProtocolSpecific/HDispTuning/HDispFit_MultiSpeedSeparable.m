%-----------------------------------------------------------------------------------------------------------------------
%-- HDispFit_MultiSpeed.m -- Fits disparity tuning curves obtained at multiple (for starters, two) speeds 
%-- with independed Gabors or with Gabors that have certain parameters tied together for different speeds.
%-- NOTE: responses are square-rooted in the error functions, to help homogenize
%-- variance (a la Cumming).  A constrained version of Levenberg-Marquardt optimization is used, as
%-- implemented using 'fmincon'.
%--	GCD, 5/23/01
%-----------------------------------------------------------------------------------------------------------------------
function HDispFit_MultiSpeedSeparable(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

% not implemented yet to select output
output = 0;

symbols = {'bo' 'r*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'b-' 'r-' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

%get the column of values of horiz. disparity in the dots_params matrix
hor_disp = data.dots_params(DOTS_HDISP,:,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (hor_disp == data.one_time_params(NULL_VALUE)) );

%get the column of speed values
speed = data.dots_params(DOTS_SPEED,:,PATCH1);
unique_speed = munique(speed(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (hor_disp == LEYE_CONTROL) | (hor_disp == REYE_CONTROL) | (hor_disp == UNCORR_CONTROL) );

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hor_disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

% Calculate spontaneous rates before looping through so can calculate DTI
null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));

if (length(unique_speed)>2) 
    disp(' more than two speeds not handled yet');
    return;
end
if (length(unique_speed)<2) 
    disp('This analysis only for runs with two speeds interleaved');
    return;
end

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'HDisp Tuning Curves fit with Gabors');

%**********************************************************************************************
%*** First, fit each data at each speed with an independent Gabor function
%**********************************************************************************************
for i=1:length(unique_speed)	
    speed_select = logical( (speed == unique_speed(i)) );
    
    plot_x = hor_disp(speed_select & ~null_trials & ~control_trials & select_trials);
    plot_y = spike_rates(speed_select & ~null_trials & ~control_trials & select_trials); 
    
    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 0);
    
    %compute the Disparity Discrimination Index
    [DDI(i), var_term] = Compute_DDI(plot_x, plot_y);
    
    %store mean rates for output
    mean_rates(i, : ) = py';
    
    p_value(i) = spk_anova(plot_y, plot_x, px);
    avg_resp(i) = mean(plot_y);      
    
    %now, compute the monoc. and uncorrelated control values
    Leye_trials = logical( (hor_disp == LEYE_CONTROL) );
    Leye_resp = spike_rates(speed_select & select_trials & Leye_trials);
    Reye_trials = logical( (hor_disp == REYE_CONTROL) );
    Reye_resp = spike_rates(speed_select & select_trials & Reye_trials);
    Uncorr_trials = logical( (hor_disp == UNCORR_CONTROL) );
    Uncorr_resp = spike_rates(speed_select & select_trials & Uncorr_trials);
    cont_rate(i).left = mean(Leye_resp);
    cont_rate(i).right = mean(Reye_resp);
    cont_rate(i).uncorr = mean(Uncorr_resp);
    
    means{i} = [px py];
    raw{i} = [plot_x' plot_y'];
    stderr{i} = perr;
    
    %plot the raw data
    subplot(2, 2, 2);
    hold on;
    PlotRawData(plot_x', plot_y', symbols{i}, max(px), Leye_resp, Leye_trials, Reye_resp, Reye_trials, Uncorr_resp, Uncorr_trials, null_rate);
    subplot(2, 2, 4);
    hold on;
    PlotRawData(plot_x', plot_y', symbols{i}, max(px), Leye_resp, Leye_trials, Reye_resp, Reye_trials, Uncorr_resp, Uncorr_trials, null_rate);
    hold off;
    
    fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
    fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
    [temp, freq(i)] = gaborfit(means{i},raw{i},fixed_param_flags,fixed_param_values);
    
    %adjust phase to be within -PI->PI
    if (temp(6) < -2*pi)
        temp(6) = temp(6) + 2*pi;
    end
    if (temp(6) > 2*pi)
        temp(6) = temp(6) - 2*pi;
    end
    pars(i,:) = temp';
    
    x_interp = (px(1)): .01 : (px(length(px)));
    y_gabor = gaborfunc(x_interp, pars(i,:));
    y_gabor(y_gabor < 0) = 0;
    
    Gabor_fit{i} = y_gabor;
    
    %compute max and min of Gabor fit
    [max_val, max_indx] = max(y_gabor);
    max{i}.rsp = max_val;
    max{i}.disp = x_interp(max_indx);
    [min_val, min_indx] = min(y_gabor);
    min{i}.rsp = min_val;
    min{i}.disp = x_interp(min_indx);
    
    % Now plot fitted curves; square the fitted values to counteract sqrt() above
    subplot(2, 2, 2);
    hold on;
    plot(x_interp, y_gabor, lines{i}, 'LineWidth', 2);
    %plot(x_interp, y_gauss, [lines{i} '-']);
    %plot(x_interp, y_sine, [lines{i} '.']);
    Gabor_SSE(i) = gaborerr(pars(i,:));
    buff = sprintf('%8.5f ', Gabor_SSE );
    title(['error = ' buff]);
    hold off;   
    
    y_fit{i} = gaborfunc(px, pars(i,:));
    y_fit{i}(y_fit{i} < 0) = 0;
    %add a column of ones to yfit to make regress happy
    y_fit{i} = [ones(length(y_fit{i}),1) y_fit{i}];
    [b, bint, r, rint, stats1(i,:)] = regress(py, y_fit{i});
    
    y_fit_raw{i} = gaborfunc(plot_x', pars(i,:));
    y_fit_raw{i}(y_fit_raw{i} < 0) = 0;
    y_fit_raw{i} = [ones(length(y_fit_raw{i}),1) y_fit_raw{i}];
    [b, bint, r, rint, stats2(i,:)] = regress(plot_y', y_fit_raw{i});
    
end
yl = YLim;
YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
XLabel('Horizontal Disparity(deg)');
YLabel('Response (spikes/sec)');

%print fit parameters to the screen
subplot(2,2,1);
axis([0 100 0 100]);    axis('off');
xpos = -20; ypos = 110;
font_size = 8;  bump_size = 9;
temp = strcat(PATH, FILE);
temp(temp == '\') = '/';
% this prevents a stupid error from appearing on the screen
line = sprintf('File: %s', temp);
text(xpos-15,ypos, line,'FontSize',font_size+1);		ypos = ypos - bump_size;
line = sprintf('Base Rate: %8.4f ', pars(:,1) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Amplitude: %8.4f ', pars(:,2));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Gauss Ctr: %8.4f ', pars(:,3));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Gauss  SD: %8.4f ', pars(:,4));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Sine Freq: %8.4f ', pars(:,5));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Sin Phase: %8.4f ', pars(:,6));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('FT Freq: %8.4f ', freq);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Means Rsq: %5.3f ', stats1(:,1) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Means P:   %6.5f ', stats1(:,3) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Raw Rsq: %5.3f ', stats2(:,1) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Raw P:   %6.5f ', stats2(:,3) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

outstr1 =  sprintf('%s %8.6f %8.6f %8.4f %8.4f %8.3f %8.3f %8.4f %8.4f %8.3f %8.3f %8.4f %8.4f %8.3f ',... 
    FILE, p_value(1), p_value(2), max{1}.disp, max{2}.disp, max{1}.rsp, max{2}.rsp, min{1}.disp, min{2}.disp, min{1}.rsp, min{2}.rsp, DDI(1), DDI(2), null_rate);
outstr2 =  sprintf('%8.3f %8.3f %8.3f %8.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.6f %8.6f %8.4f %8.4f %8.6f %8.6f ',... 
    pars(1,1), pars(2,1), pars(1,2), pars(2,2), pars(1,3), pars(2,3), pars(1,4), pars(2,4), pars(1,5), pars(2,5), pars(1,6), pars(2,6), stats1(1,1), stats1(2,1), stats1(1,3), stats1(2,3), stats2(1,1), stats2(2,1), stats2(1,3), stats2(2,3) );

%also write out data in form suitable for plotting tuning curve with Origin.
i = size(PATH,2) - 1;
while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
    i = i - 1;
end   
PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
i = size(FILE,2) - 1;
while FILE(i) ~='.'
    i = i - 1;
end
FILEOUT2 = [FILE(1:i) 'hdsp_gabor_multispeed'];
fileid = [PATHOUT FILEOUT2];
proffid = eval(['fopen(fileid, ''w'')']);
fprintf(proffid,'InDsp\tStatFit\tMovFit\tHDisp\tStatRsp\tStatErr\tMovResp\tMovErr\tHDisp2\tStatLRU\tStatLbl\tMovLRU\tMovLbl\tHDisp3\tSpon\r\n');
N_interp_disp = length(x_interp);
N_disp = length(means{1}(:,1));
for i=1:N_interp_disp
    fprintf(proffid,'%6.3f\t%6.3f\t%6.3f\t', x_interp(i), Gabor_fit{1}(i), Gabor_fit{2}(i) );
    if (i <= N_disp)
        fprintf(proffid,'%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t', means{1}(i,1), means{1}(i,2), stderr{1}(i), means{2}(i,2), stderr{2}(i) );
    else
        fprintf(proffid,'\t\t\t\t\t');
    end
    if (i == 1)
        fprintf(proffid,'%6.2f\t%6.2f\t%3s\t%6.2f\t%3s\t%6.2f\t%6.2f\r\n',1.15*means{1}(N_disp,1),cont_rate(1).left,'L',cont_rate(2).left,'L',means{1}(1,1),null_rate);
    elseif (i==2)		
        fprintf(proffid,'%6.2f\t%6.2f\t%3s\t%6.2f\t%3s\t%6.2f\t%6.2f\r\n',1.15*means{1}(N_disp,1),cont_rate(1).right,'R',cont_rate(2).right,'R',means{1}(N_disp,1),null_rate);
    elseif (i==3)		
        fprintf(proffid,'%6.2f\t%6.2f\t%3s\t%6.2f\t%3s\r\n',1.15*means{1}(N_disp,1),cont_rate(1).uncorr,'U',cont_rate(2).uncorr,'U');
    else
        fprintf(proffid,'\r\n');
    end
end
fclose(proffid);


%**********************************************************************************************
%*** Now, fit data at both speeds with a pair of Gabors having only independent amplitude and base rate
%**********************************************************************************************
[pars] = separable_gaborfit(means,raw);

subplot(2, 2, 4);
hold on;
pars1 = pars(1:6); %parameters for the 1st Gabor
gabor1 = gaborfunc_alt(x_interp, pars1);
gabor1(gabor1 < 0) = 0;
plot(x_interp, gabor1, lines{1}, 'LineWidth', 2);

hold on;
pars2 = pars(1:6); %parameters for the 1st Gabor
pars2(1) = pars(7);
gabor2 = gaborfunc_alt(x_interp, pars2);
gabor2(gabor2 < 0) = 0;
plot(x_interp, gabor2, lines{2}, 'LineWidth', 2);

paired_err = separable_gaborerr(pars);
buff = sprintf('Paired error = %8.5f ', paired_err );
title(buff);
hold off;   

%now compute an overall R^2 for the fit with same_shape Gabors
gabor1_coarse = gaborfunc_alt(means{1}(:,1), pars1);
gabor1_coarse(gabor1_coarse < 0) = 0;
gabor1_coarse = [ones(length(gabor1_coarse),1) gabor1_coarse];
[b1, bint1, r1, rint1, stats_constr1] = regress(means{1}(:,2), gabor1_coarse);
gabor2_coarse = gaborfunc_alt(means{2}(:,1), pars2);
gabor2_coarse(gabor2_coarse < 0) = 0;
gabor2_coarse = [ones(length(gabor2_coarse),1) gabor2_coarse];
[b2, bint2, r2, rint2, stats_constr2] = regress(means{2}(:,2), gabor2_coarse);

%do Sequential F-test
independent_err = sum(Gabor_SSE);
indep_params = 12;  paired_params = 7;
Npts = length(spike_rates(~null_trials & ~control_trials & select_trials));
Fseq = ( (paired_err - independent_err )/(indep_params-paired_params) ) / ( independent_err/(Npts - indep_params) );
Pseq = 1 - fcdf(Fseq, (indep_params-paired_params), (Npts-indep_params) );

%print fit parameters to the screen
subplot(2,2,3);
axis([0 100 0 100]);    axis('off');
xpos = -20; ypos = 110;
font_size = 8;  bump_size = 9;
line = sprintf('Scale 1: %8.4f ', pars(1) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Ampl/BRate: %8.4f ', pars(2));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Gauss Ctr: %8.4f ', pars(3));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Gauss  SD: %8.4f ', pars(4));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Sine Freq: %8.4f ', pars(5));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Sin Phase: %8.4f ', pars(6));
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Scale 2: %8.4f ', pars(7) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

ypos = ypos - bump_size;
line = sprintf('Sequential F-test:');
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('F = %7.3f, P = %7.5f', Fseq, Pseq);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;


outstr3 = sprintf('%8.4f %8.6f %8.6f %8.6f', Fseq, Pseq, stats_constr1(1), stats_constr2(1));

outfile = [BASE_PATH 'ProtocolSpecific\HDispTuning\MultiSpeed_GaborFit_Separable.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, 'FILE          AnovP1   AnovP2   Dmax1   Dmax2    Rmax1   Rmax2    Dmin1   Dmin2    Rmin1   Rmin2    DDI1    DDI2     Spont    BRate1  BRate2   Amp1    Amp2     GsCtr1  GsCtr2   GsSD1   GsSD2    Freq1   Freq2    Phs1     Phs2     Rsqmean1  Rsqmean2   Pmean1  Pmean2   Rsqraw1   Rsqraw2    Praw1    Praw2     Fseq     Pseq     RsqCon1   RsqCon2 ');
    fprintf(fid, '\r\n');
end
fprintf(fid, '%s', [outstr1 outstr2 outstr3]);
fprintf(fid, '\r\n');
fclose(fid);

%now, print out some useful information in the upper subplot
%PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

return;

%----------------------------------------
function    PlotRawData(x, y, symb1, max_x, L_resp, L_trials, R_resp, R_trials, U_resp, U_trials, null_resp)

plot(x, y, symb1);

hold on;
errorbar(max_x*1.15, mean(L_resp), std(L_resp)/sqrt(sum(L_trials)), std(L_resp)/sqrt(sum(L_trials)), symb1);
text(max_x*1.3, mean(L_resp), 'L');
hold on;
errorbar(max_x*1.15, mean(R_resp), std(R_resp)/sqrt(sum(R_trials)), std(R_resp)/sqrt(sum(R_trials)), symb1);
text(max_x*1.3, mean(R_resp), 'R');
hold on;
errorbar(max_x*1.15, mean(U_resp), std(U_resp)/sqrt(sum(U_trials)), std(U_resp)/sqrt(sum(U_trials)), symb1);
text(max_x*1.3, mean(U_resp), 'U');

%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
null_x = [min(x) max(x)];
null_y = [null_resp null_resp];
hold on;
plot(null_x, null_y, 'k--');
hold off;

return;