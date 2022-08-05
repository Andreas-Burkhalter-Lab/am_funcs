function roi_traceviewer(basename, intensityfilename, options)
%roi_traceviewer takes the roi intensity matrix calculated by
%calc_stack_roi_inten and returns a set of panels with useful plots. As the
%user hits any key the panels refresh to show information from the next
%sequential ROI. The raw data, a bessel filtered version of the data, and
%trial overlays for each stimulus (plotted using trialcolor for
%color picking)as well as the mean response for all trials.
% Syntax:
%   roi_traceviewer(basename, intensityfilename)
%   roi_traceviewer(basename, intensityfilename, options)
% where
%   basename is the original basic filename from the imagine experimental data (also the header name), 
%             BTW:  the header file from the original data should automatically be saved
%             as a .mat or .imagine file into the parent directory for this job by
%             the previous processing function (registration, etc.). If not, you should do this.
%             NOTE: if the previous processing function modified the shape
%             of the data (e.g. dropped the first and last stacks) the
%             copied header needs to be edited to reflect this (you may need to check that
%             this is so).  
%   intensityfilename is the name of the saved intensity matrix file created by
%             calc_stack_roi_inten 
%   options is a structure indicating some of the extra behavior you would
%             or would not like roi_traceviewer to exhibit.
%        .timewindowpretrans can be set to change the pre-stimulus period
%             used for the trial snippeting time window.
%             Default is 4 (stacks).
%        .timewindowposttrans can be set to change the after stimulus-transition period
%             used for the trial snippeting time window.
%             Default is 15 (stacks).
%        .backgrndmeanlength can be set to change the pre-stimulus background period
%             used for normalization.
%             Default is 3 (stacks).
%        .n is the number of poles to used for the bessel noise reduction
%             filter. Default is 3.
%        .Wn Wn is a two-element vector, Wn = [w1 w2], butter returns an order 2*n 
%             digital bandpass filter with passband w1 < w < w2. Default is [0.02 0.5].
%        .filt_type indicates the type of filter you want to use. Options
%             are: 'high', 'low', 'stop'.  If you say nothing, the automatic
%             default is bandpass.
%NOTE: for the moment we only have this set up to handle the short movies
%that have only 4 stimulus types. Mods will be required for longer, bigger,
%better experimental data.
% 
%Copywrite 2006 by Diwakar Turaga & Terrence Holekamp

if nargin < 3
    options.timewindowpretrans = 4;
    options.timewindowposttrans = 15;
    options.backgrndmeanlength = 3;
    options.n = 3;
    options.Wn = [0.1 0.5];
end

if(~isfield(options, 'timewindowpretrans'))
    options.timewindowpretrans = 4;
end
if(~isfield(options, 'timewindowposttrans'))
    options.timewindowposttrans = 15;
end
if(~isfield(options, 'backgrndmeanlength'))
    options.backgrndmeanlength = 3;
end
if(~isfield(options, 'n'))
    options.n = 3;
end
if(~isfield(options, 'Wn'))
    options.Wn = [0.1 0.5];
end

timewindowpretrans = options.timewindowpretrans;
timewindowposttrans = options.timewindowposttrans;
bkglng = options.backgrndmeanlength;
n = options.n;
Wn = options.Wn;

%[b,a] = besself(n,Wo);
if(isfield(options, 'filt_type'))
    [b,a] = butter(n,Wn,options.filt_type);
else
    [b,a] = butter(n,Wn);
end

if ~exist('intensities','var')
    load (intensityfilename,'-mat')
end

smm = stackmm(basename);
header = smm.header;
stimuli = header.stim_lookup;

size_inten = size(intensities);
number_of_ROI  = size_inten(1);
number_of_stacks = size_inten(2);
number_of_trials = max(unique(header.trial_lookup(~isnan(header.trial_lookup))));

intensities(:,1) = intensities(:,2);
intensities = intensities';
%intensities_filt = filter(b,a,intensities);
intensities_filt = filtfilt(b,a,intensities);
intensities_filt = intensities_filt';
intensities = intensities';

% for i = 1:number_of_ROI
%     inten_avg = nanmean(intensities_filt(i,:));
%     inten_normal(i,:) = (intensities_filt(i,:) - inten_avg)/inten_avg;
% end

for i = 1:number_of_ROI
    inten_avg_unfilt = nanmean(intensities(i,:));
    inten_normal_unfilt(i,:) = (intensities(i,:) - inten_avg_unfilt)/inten_avg_unfilt;
end

if bkglng > timewindowpretrans
    disp('you may be asking for too long a background averaging period')
end


figure;

stimulus_1 = NaN(number_of_trials,number_of_ROI, timewindowpretrans+timewindowposttrans+1);
stimulus_2 = NaN(number_of_trials,number_of_ROI, timewindowpretrans+timewindowposttrans+1);
stimulus_3 = NaN(number_of_trials,number_of_ROI, timewindowpretrans+timewindowposttrans+1);
stimulus_4 = NaN(number_of_trials,number_of_ROI, timewindowpretrans+timewindowposttrans+1);

i = 1;
trial_counter1 = 1; trial_counter2 = 1; trial_counter3 = 1; trial_counter4 = 1;
while i <=number_of_stacks-1
    
    if (stimuli(i)==0 && stimuli(i+1)~=0)
            
            if stimuli(i+1) == 1
                for k = 1:number_of_ROI
                    inten_backgrnd_avg1 = nanmean(intensities(k,i-bkglng:i));
                    if i+timewindowposttrans <= size(intensities,2)
                        inten_snip1 = intensities(k,i-timewindowpretrans:i+timewindowposttrans);
                    else
                        inten_snip1 = intensities(k,i-timewindowpretrans:end);
                    end
                    inten_normal1 = (inten_snip1 - inten_backgrnd_avg1)/inten_backgrnd_avg1;
                    stimulus_1(trial_counter1,k,1:length(inten_normal1)) = inten_normal1;
                end
                i = i+1;
                trial_counter1 = trial_counter1+1;
            end   
            
            if stimuli(i+1) == 2
                for k = 1:number_of_ROI
                    inten_backgrnd_avg2 = nanmean(intensities(k,i-bkglng:i));
                    if i+timewindowposttrans <= size(intensities,2)
                        inten_snip2 = intensities(k,i-timewindowpretrans:i+timewindowposttrans);
                    else
                        inten_snip2 = intensities(k,i-timewindowpretrans:end);
                    end
                    inten_normal2 = (inten_snip2 - inten_backgrnd_avg2)/inten_backgrnd_avg2;
                    stimulus_2(trial_counter2,k,1:length(inten_normal2)) = inten_normal2;
                end
                i = i+1;
                trial_counter2 = trial_counter2+1;
            end 
            
            if stimuli(i+1) == 3
                for k = 1:number_of_ROI
                    inten_backgrnd_avg3 = nanmean(intensities(k,i-bkglng:i));
                    if i+timewindowposttrans <= size(intensities,2)
                        inten_snip3 = intensities(k,i-timewindowpretrans:i+timewindowposttrans);
                    else
                        inten_snip3 = intensities(k,i-timewindowpretrans:end);
                    end
                    inten_normal3 = (inten_snip3 - inten_backgrnd_avg3)/inten_backgrnd_avg3;
                    stimulus_3(trial_counter3,k,1:length(inten_normal3)) = inten_normal3;
                end
                i = i+1;
                trial_counter3 = trial_counter3+1;
            end 
            
            if stimuli(i+1) == 4
                for k = 1:number_of_ROI
                    inten_backgrnd_avg4 = nanmean(intensities(k,i-bkglng:i));
                    if i+timewindowposttrans <= size(intensities,2)
                        inten_snip4 = intensities(k,i-timewindowpretrans:i+timewindowposttrans);
                    else
                        inten_snip4 = intensities(k,i-timewindowpretrans:end);
                    end
                    inten_normal4 = (inten_snip4 - inten_backgrnd_avg4)/inten_backgrnd_avg4;
                    stimulus_4(trial_counter4,k,1:length(inten_normal4)) = inten_normal4;
                end
                i = i+1;
                trial_counter4 = trial_counter4+1;
            end        
            
    else
        i=i+1;
    end
end

for i = 1:number_of_ROI
    
    meanstim1 = squeeze(nanmean(stimulus_1(:,i,:),1));
    meanstim2 = squeeze(nanmean(stimulus_2(:,i,:),1));
    meanstim3 = squeeze(nanmean(stimulus_3(:,i,:),1));
    meanstim4 = squeeze(nanmean(stimulus_4(:,i,:),1));
    
    subplot(2,3,1), plot(intensities_filt(i,:)), title(['ROI # ' num2str(i)]);
    subplot(2,3,4), plot(inten_normal_unfilt(i,:)), title('unfilt norm data');

    stimulus_1_trialcount = size(stimulus_1,1);
    trialmap1 = trialcolor(stimulus_1_trialcount);
    for j = 1:stimulus_1_trialcount    
        subplot(2,3,2), plot(squeeze(stimulus_1(j, i, :)),'color',trialmap1{j}), axis tight, hold on,title('stimulus 1');
    end
    plot(meanstim1,'color',[0 0 0],'LineWidth',2)
    hold off
    xlimits = get(gca,'xlim');
    ylimits = get(gca,'ylim');
    
    stimulus_2_trialcount = size(stimulus_2,1);
    trialmap2 = trialcolor(stimulus_2_trialcount);
    for j = 1:stimulus_2_trialcount
        subplot(2,3,3), plot(squeeze(stimulus_2(j, i, :)),'color',trialmap2{j}), axis tight, hold on,title('stimulus 2');
    end
    plot(meanstim2,'color',[0 0 0],'LineWidth',2)
    set(gca, 'xlim', xlimits);
    set(gca, 'ylim', ylimits);
    hold off
    
    stimulus_3_trialcount = size(stimulus_3,1);
    trialmap3 = trialcolor(stimulus_3_trialcount);
    for j = 1:stimulus_3_trialcount
        subplot(2,3,5), plot(squeeze(stimulus_3(j, i, :)),'color',trialmap3{j}), axis tight, hold on,title('stimulus 3');
    end
    plot(meanstim3,'color',[0 0 0],'LineWidth',2)
    set(gca, 'xlim', xlimits);
    set(gca, 'ylim', ylimits);
    hold off

    stimulus_4_trialcount = size(stimulus_4,1);
    trialmap4 = trialcolor(stimulus_4_trialcount);
    for j = 1:stimulus_4_trialcount
        subplot(2,3,6), plot(squeeze(stimulus_4(j, i, :)),'color',trialmap4{j}), axis tight, hold on, title('stimulus 4');
    end
    plot(meanstim4,'color',[0 0 0],'LineWidth',2)
    set(gca, 'xlim', xlimits);
    set(gca, 'ylim', ylimits);
    hold off

    pause;
    clf;

end


%
% for i = 1:number_of_ROI; 
%     clear inten_normal;
%     inten_total = 0;
%     
%     for j = 1:Number_to_avg
%         inten_total = inten_total + inten(i, j);
%     end
%     
%     inten_avg = inten_total/Number_to_avg;
%     
%     inten_normal(i,:) = (inten(i,:) - inten_avg)./inten_avg;
%     
%     inten_normal(i,1) = inten_normal(i,2);
%     subplot(2,1,1), plot(inten(i,:))
%     subplot(2,1,2), plot(stim) 
%     title(num2str(i));
%     pause;
%     
% 
% end