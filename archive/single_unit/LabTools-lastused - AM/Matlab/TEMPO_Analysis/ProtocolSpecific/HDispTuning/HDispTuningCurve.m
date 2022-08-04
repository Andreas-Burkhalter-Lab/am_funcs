%-----------------------------------------------------------------------------------------------------------------------
%-- HDispTuningCurve.m -- Plots a horizontal disparity tuning curve, possibly for multiple speeds, and including
%--	monoc and uncorrelated control conditions.
%--	GCD, 1/26/00
%-----------------------------------------------------------------------------------------------------------------------
function HDispTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

%get the column of values of horiz. disparity in the dots_params matrix
hor_disp = data.dots_params(DOTS_HDISP,:,PATCH1);
%hor_disp'

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (hor_disp == data.one_time_params(NULL_VALUE)) );

%get the column of speed values
speed = data.dots_params(DOTS_SPEED,:,PATCH1);
unique_speed = munique(speed(~null_trials)')

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%compute the horizontal vergence angle for each trial
h_verg = data.eye_positions(LEYE_H, :) - data.eye_positions(REYE_H, :);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (hor_disp == LEYE_CONTROL) | (hor_disp == REYE_CONTROL) | (hor_disp == UNCORR_CONTROL) );

unique_hdisp = munique(hor_disp(~null_trials & ~control_trials)')

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hor_disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%calculate the average intra-trial variance in vergence angle, added by GCD 9/26/02
verg_trial_std = [];
eye_traces = ComputeAllEyeTraces(data, length(hor_disp), StartCode, StopCode, 0, 0);
for i = 1:length(hor_disp)  %loop through all trials
    h_verg_trace = eye_traces{LEYE_H,i} - eye_traces{REYE_H,i};
    verg_trial_std(i) = std(h_verg_trace);
end
avg_intratrial_verg_std = mean(verg_trial_std);
avg_intertrial_verg_std = std(h_verg);

% Calculate spontaneous rates before looping through so can calculate DTI
null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
null_rate = mean(null_resp);

figure(2);
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [450 50 500 573], 'Name', 'Horizontal Disparity Tuning Curve');
subplot(2, 1, 2);

%for i=1:length(unique_speed)
%    disptune{i} = sprintf('%s%d%s', 'disptuning', i, '.mat');
%    plot_val{i} = sprintf('%s%d%s', 'plot_values', i, '.mat');
%end

for i=1:length(unique_speed)	%for each different speed value, plot a separate disparity tuning curve
    speed_select = logical( (speed == unique_speed(i)) );
    
    plot_x = hor_disp(speed_select & ~null_trials & ~control_trials & select_trials);
    plot_y = spike_rates(speed_select & ~null_trials & ~control_trials & select_trials); 
    verg_select = h_verg(speed_select & ~null_trials & ~control_trials & select_trials);
    
    %compute the Disparity Discrimination Index
    [DDI(i), var_term(i)] = Compute_DDI(plot_x, plot_y);
    
    %compute the tuning curve Asymmetry Index
    ASI(i) = Compute_ASI(plot_x, plot_y);
    
    %Do an ANCOVA to test signif of disparity tuning in presence of vergence
    [H,ATAB,CTAB,STATS] = aoctool(verg_select, plot_y, plot_x);
    AncPdisp(i) = ATAB{2,6};
    AncPverg(i) = ATAB{3,6};
    AncPinter(i) = ATAB{4,6};
    close(3); close(4); close(5);

    verg_sd(i) = std(verg_select);
    
    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    figure(2);
    hold on;
    [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 1);
        
    %Compute DTI from spline fit
    DTI(i) = 1 - (pmin(i).y - null_rate)/(pmax(i).y - null_rate);

    %Calculate modulation index using sqrt raw responses and subtracting spontaneous
    DMI(i) = Compute_ModIndex(plot_x, plot_y, null_resp)
    
    %Calculate ROC from max and min responses TU 12/19/01
    [max_value, max_index] = max(py);
    [min_value, min_index] = min(py);
    pref_dist = plot_y(logical(plot_x == px(max_index)));
    null_dist = plot_y(logical(plot_x == px(min_index)));
    ROC_value(i) = rocN(pref_dist, null_dist, 100);
    
    %[px py perr]
    
    %store mean rates for output
    mean_rates(i, : ) = py';
    
    hold off;
    
    hold on;
    errorbar(px, py, perr, perr, symbols{i});
    hold off;
    
    % do the ANOVA on sqrt(firing rate) a la Prince et al.
    [p_value(i), MSgroup(i), MSerror(i)] = spk_anova_F(sqrt(plot_y), plot_x, px);
    avg_resp(i) = mean(plot_y);      
    
    %now, compute and plot the monoc. and uncorrelated control values
    Leye_trials = logical( (hor_disp == LEYE_CONTROL) );
    Leye_resp = spike_rates(speed_select & select_trials & Leye_trials);
    hold on;
    errorbar(max(px)*1.07, mean(Leye_resp), std(Leye_resp)/sqrt(sum(Leye_trials)), std(Leye_resp)/sqrt(sum(Leye_trials)), symbols{i});
    text(max(px)*1.12, mean(Leye_resp), 'L');
    
    Reye_trials = logical( (hor_disp == REYE_CONTROL) );
    Reye_resp = spike_rates(speed_select & select_trials & Reye_trials);
    hold on;
    errorbar(max(px)*1.07, mean(Reye_resp), std(Reye_resp)/sqrt(sum(Reye_trials)), std(Reye_resp)/sqrt(sum(Reye_trials)), symbols{i});
    text(max(px)*1.12, mean(Reye_resp), 'R');
    
    Uncorr_trials = logical( (hor_disp == UNCORR_CONTROL) );
    Uncorr_resp = spike_rates(speed_select & select_trials & Uncorr_trials);
    hold on;
    errorbar(max(px)*1.07, mean(Uncorr_resp), std(Uncorr_resp)/sqrt(sum(Uncorr_trials)), std(Uncorr_resp)/sqrt(sum(Uncorr_trials)), symbols{i});
    text(max(px)*1.12, mean(Uncorr_resp), 'U');
    
    %added so that the grad model can read in the disp tuning data
    %save_uncorr_resp = mean(Uncorr_resp);
    %save(disptune{i}, 'px', 'py', 'perr', 'save_uncorr_resp');
    %save disptune px py perr save_uncorr_resp;
    
    cont_rate(i).left = mean(Leye_resp);
    cont_rate(i).right = mean(Reye_resp);
    cont_rate(i).uncorr = mean(Uncorr_resp);
end

corr_coef = corrcoef(mean_rates');
if size(corr_coef,1) > 1
    corr_coef = corr_coef(2);
else
    corr_coef = NaN;
end  
%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
null_x = [min(px) max(px)];
null_y = [null_rate null_rate];
hold on;
plot(null_x, null_y, 'k--');
hold off;

yl = YLim;
YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
XLabel('Horizontal Disparity(deg)');
YLabel('Response (spikes/sec)');

%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%now print out useful values for H.Disp specific 
% pmax, pmin, py, 
PrintHDispData(p_value, avg_resp, pmax, pmin, px, null_rate, cont_rate, unique_speed, PATH, FILE, DDI, DTI, corr_coef, ASI);

output = 1;
if (output == 1)
    i = size(PATH,2) - 1;
    while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
        i = i - 1;
    end   
    PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
    i = size(FILE,2) - 1;
    while FILE(i) ~='.'
        i = i - 1;
    end
    FILEOUT = [FILE(1:i) 'hdsp'];
    
    fileid = [PATHOUT FILEOUT];
%    fwriteid = eval(['fopen(fileid, ''w'')']);
    % header for data fields
    %fprintf(fwriteid, '%%/speed /	ave resp	/spont resp	/max resp	/hdisp at max/	min resp/	hdisp at min/	Curve ANOVA/DTI/	DDI / CorrCoef	/Disp -1.6	Disp -1.2	Disp -0.8	Disp -0.4	Disp 0	Disp +0.4	Disp +0.8	Disp +1.2	Disp +1.6	L Control	R Control	U Control	Pref Dir	Pref Speed	Mapped Pref H Disp	RF X-Ctr	RF Y-Ctr	RF Diam\n');
    
%    for speed = 1:length(unique_speed)
%        fprintf(fwriteid, '%3.1f	%5.2f	%5.2f	%6.2f	%5.2f	%5.2f	%5.2f	%3.2f %3.2f %3.2f %3.2f %5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f\n', unique_speed(speed), avg_resp(speed), null_rate, pmax(speed).y, pmax(speed).x, pmin(speed).y, pmin(speed).x, p_value(speed), DTI(speed), DDI(speed), corr_coef, mean_rates(speed,1), mean_rates(speed,2), mean_rates(speed,3), mean_rates(speed,4), mean_rates(speed,5), mean_rates(speed,6), mean_rates(speed,7), mean_rates(speed,8), mean_rates(speed,9), cont_rate(speed).left, cont_rate(speed).right, cont_rate(speed).uncorr, data.one_time_params(PREFERRED_DIRECTION), data.one_time_params(PREFERRED_SPEED), data.one_time_params(PREFERRED_HDISP), data.one_time_params(RF_XCTR), data.one_time_params(RF_YCTR), data.one_time_params(RF_DIAMETER));	
%    end    
%    fclose(fwriteid);
    
    %----------------------------------------------------------------------------
    %also write out data in form suitable for plotting tuning curve with Origin.
    FILEOUT2 = [FILE(1:i) 'hdsp_curv'];
    fileid = [PATHOUT FILEOUT2];
    proffid = eval(['fopen(fileid, ''w'')']);
    
    fprintf(proffid,'HDisp\tAvgResp\tStdErr\tHDisp2\tMonoc\tMonErr\tLabel\tHDisp3\tSpon\n');
    for i=1:length(px)
        fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t', px(i), py(i), perr(i));
        if (i == 1)
            fprintf(proffid,'%6.2f\t%6.2f\t%6.2f\t%6.1f\t%6.2f\t%6.2f\n',1.15*px(length(px)),mean(Leye_resp),std(Leye_resp)/sqrt(sum(Leye_trials)),-99,null_x(i),null_y(i));
        elseif (i==2)		
            fprintf(proffid,'%6.2f\t%6.2f\t%6.2f\t%6.1f\t%6.2f\t%6.2f\n',1.15*px(length(px)),mean(Reye_resp),std(Reye_resp)/sqrt(sum(Reye_trials)),99,null_x(i),null_y(i));
        elseif (i==3)		
            fprintf(proffid,'%6.2f\t%6.2f\t%6.2f\t%6.1f\t%6.2f\t%6.2f\n',1.15*px(length(px)),mean(Uncorr_resp),std(Uncorr_resp)/sqrt(sum(Uncorr_trials)),98,null_x(1),null_y(1));
        else
            %				fprintf(proffid,'*\t*\t*\t*\t*\t*\n');
            fprintf(proffid,'\t\t\t\t\t\n');
        end
    end
    
    fclose(proffid);
    %----------------------------------------------------------------------------
    
    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimulus speed for each neuron.
    %outfile = ['Z:\Users\jerry\GradAnalysis\HDispTuningSummary_8.27.02.dat'];
    outfile = [BASE_PATH 'ProtocolSpecific\HDispTuning\HDispTuningSummary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Spd\t MaxDsp\t MaxRsp\t MinDsp\t MinRsp\t AvgRsp\t Spont\t DspMI\t Leye\t Reye\t Uncorr\t DTI\t AnovaP\t DDI\t VarTerm\t ROC\t MSgroup\t MSerror\t AncovaPDisp\t AncovaPverg\t AncovaPint\t VergSD\t LgDsp\t SmDsp\t VergSDIntra\t VergSDInter\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    for j = 1:length(unique_speed)
        buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %5.2f\t %6.3f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %10.8f\t %6.3f\t %6.3f\t %5.3f\t %10.3f\t %10.3f\t %10.8f\t %10.8f\t %10.8f\t %6.4f\t %6.4f\t %6.4f\t %10.8f\t %10.8f\t', ...
            FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
            unique_speed(j), pmax(j).x, pmax(j).y, pmin(j).x, pmin(j).y, avg_resp(j), null_rate, DMI(j), cont_rate(j).left, cont_rate(j).right, cont_rate(j).uncorr, DTI(j), p_value(j), DDI(j), var_term(j), ROC_value(j), MSgroup(j), MSerror(j), AncPdisp(j), AncPverg(j), AncPinter(j), verg_sd(j), max(unique_hdisp), min(unique_hdisp), avg_intratrial_verg_std, avg_intertrial_verg_std );
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
    end
    fclose(fid);
    %------------------------------------------------------------------------
    
    %print out data for gradient use.  This is now redundant.  GCD 8/8/01    
    grad_print = 0;
    if (grad_print == 1)
        %pathsize = size(PATH,2) - 1;
        %while PATH(pathsize) ~='\'	%Analysis directory is one branch below Raw Data Dir
        %    pathsize = pathsize - 1;
        %end   
        PATHOUT = 'Z:\Data\Tempo\GradAnalysis\';
        
        for speed = 1:length(unique_speed)
            line = sprintf('%s %1.12f	%1.12f', FILE, DTI(speed), DDI(speed));	
        end    
        
        outfile = [PATHOUT 'HDisp_metrics.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE         DTI                DDI');
            fprintf(fid, '\r\n');
        end
        fprintf(fid, '%s', [line]);
        fprintf(fid, '\r\n');
        fclose(fid);
    end
    
    
end  %if (output == 1)

return;