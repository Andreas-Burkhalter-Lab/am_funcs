function roi_traceviewer_2(basename, intensityfilename, options)
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
%        .n is the number of poles to used for the butterworth noise reduction
%             filter. Default is 5.
%        .Wn is a one or two-element vector, Wn = [w1 w2], butter returns an order 2*n 
%             digital bandpass filter with passband w1 < w < w2. Default is lowpass [0.45].
%        .filt_type indicates the type of filter you want to use. Options
%             are: 'high', 'low', 'stop'.  Default is lowpass.
%        .filterpaneldata indicates whether you want the filtered or
%             unfiltered data presented in the trial overlay panels.
%             Default is 1(true) for presenting the filtered data.
%
%Copywrite 2006 by Diwakar Turaga & Terrence Holekamp

if nargin < 3
    options.timewindowpretrans = 4;
    options.timewindowposttrans = 15;
    options.backgrndmeanlength = 3;
    options.n = 5;
    options.Wn = [0.45];
    options.paneldatatype = 1;
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
    options.n = 5;
end
if(~isfield(options, 'Wn'))
    options.Wn = [0.45];
end
if(~isfield(options, 'filterpaneldata'))
    options.filterpaneldata = 1;
end

timewindowpretrans = options.timewindowpretrans;
timewindowposttrans = options.timewindowposttrans;
bkglng = options.backgrndmeanlength;
n = options.n;
Wn = options.Wn;
filtpanel = options.filterpaneldata;

if(isfield(options, 'filt_type'))
    [b,a] = butter(n,Wn,options.filt_type);
else
    [b,a] = butter(n,Wn,'low');
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
% number_of_trials = max(unique(header.trial_lookup(~isnan(header.trial_lookup))));
number_of_diff_stimuli = max(unique(stimuli));

number_of_trials = 0; % this is so that if we dont register the whole stacks, we can find out exactly how many of the trials had been registered
for i = 1:number_of_stacks
    if stimuli(i) == 1 && stimuli(i-1) == 0
        number_of_trials = number_of_trials + 1;
    end
end

% filter the data
intensities(:,1) = intensities(:,2);
intensities_filt = filtfilt(b,a,intensities');
intensities_filt = intensities_filt';

header = imreadheader(basename);
ranges = segment(header.stim_lookup, 0);
ranges = reshape(ranges,2,[]);
ranges = ranges';

for i = 1:number_of_ROI
    inten_avg_unfilt = nanmean(intensities(i,:));
    inten_normal_unfilt(i,:) = (intensities(i,:) - inten_avg_unfilt)/inten_avg_unfilt;
end

if bkglng > timewindowpretrans
    disp('you may be asking for too long a background averaging period')
end

if filtpanel
    data = intensities_filt;
else
    data = intensities;
end

stimulus = NaN(number_of_trials,number_of_ROI, timewindowpretrans+timewindowposttrans+1, number_of_diff_stimuli);

figure;
for i = 1:number_of_diff_stimuli
    trial_counter(i) = 1;
end

i = 1;

while i <=number_of_stacks-1
    
    if (stimuli(i)==0 && stimuli(i+1)~=0)
            
        for j = 1:number_of_diff_stimuli
           
            if stimuli(i+1) == j
                for k = 1:number_of_ROI
                    inten_backgrnd_avg = nanmean(data(k,i-bkglng:i));
                    if i+timewindowposttrans <= size(data,2)
                        inten_snip = data(k,i-timewindowpretrans:i+timewindowposttrans);
                    else
                         inten_snip = data(k,i-timewindowpretrans:end);
                    end
                    inten_normal = (inten_snip - inten_backgrnd_avg)/inten_backgrnd_avg;
                    stimulus(trial_counter(j),k,1:length(inten_normal),j) = inten_normal;
                end
                i = i+1;
                trial_counter(j) = trial_counter(j)+1;
            end 
            
        end
                        
            
    else
        i=i+1;
    end
end

total_plots = number_of_diff_stimuli + 4; % raw and filt stimuli take up 2 each

num_rows = ceil(total_plots/4);

ylimits = nan(2,number_of_diff_stimuli-1);
hax = nan(1,number_of_diff_stimuli);

for i = 1:number_of_ROI
    
    for j = 1: number_of_diff_stimuli
        meanstim(j,:) = (squeeze(nanmean(stimulus(:,i,:,j),1)))';
    end


    subplot(num_rows,4,[1 2]), plot(intensities_filt(i,:),'color',[0 0 0]), title('Lowpass Filtered Data');
    subplot(num_rows,4,[3 4]), plot(inten_normal_unfilt(i,:)), title('Unfiltered Normalized Data');

    for j = 1:number_of_diff_stimuli
        r = diff(find(header.stim_lookup ==1));
        r(r == 1) = 0;
        s = find(r);
        stimrange = s(1);
        
        trialmap = trialcolor(number_of_trials);
        for k = 1:number_of_trials    
            hax(j) = subplot(num_rows,4,j+4); 
            hold on
            plot((squeeze(stimulus(k, i, :, j)))','color',trialmap{k});
            axis tight
            set(gca,'TickDir','out','box','off');
            title(header.stim_labels{j});    
        end
        if j > 1
            ylimits(:,j-1)  = get(gca,'ylim');
        end
        
        plot(meanstim(j,:),'color',[0 0 0],'LineWidth',2.5)
        yaxrange = get(gca,'Ylim');
        xaxrange = get(gca,'Xlim');
        tracehand = get(gca,'Children');
        hrect(j) = rectangle('Position',[timewindowpretrans+1, yaxrange(1), stimrange,...
            yaxrange(2)-yaxrange(1)],'EdgeColor','none','FaceColor',[0.97 0.97 0.97]);
        set(gca,'Xlim',xaxrange,'Ylim',yaxrange);
        set(gca,'Children',[tracehand;hrect(j)']);
        hold off
%         legend('toggle')
    end
    yl = abs(ylimits);
    [yuk,Ind] = max(yl,[],2);
    ylimi = [ylimits(1,Ind(1)) ylimits(2,Ind(2))];
    set(hax(2:number_of_diff_stimuli),'YLim',ylimi)
    for j = 2:number_of_diff_stimuli
        yaxrange = get(hax(j),'Ylim');
        posishun = get(hrect(j),'position');
        set(hrect(j),'position',[posishun(1) yaxrange(1) posishun(3) yaxrange(2)-yaxrange(1)]);
    end
    suptitle(['ROI # ' num2str(i)]);
    pause;
    clf;

end


