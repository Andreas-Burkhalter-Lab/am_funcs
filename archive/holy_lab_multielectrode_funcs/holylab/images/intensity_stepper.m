function mfig = intensity_stepper(inten,roi,varargin)
% intensity_stepper(intenfile,roifile) opens a keystepper-controlled window of your ROIs
% TODO: help

% Copyright 2010-2011 Julian P Meeks (Timothy Holy Laboratory)

% TODO: logically parse inputs
if isstruct(inten)
	intens = inten; clear inten;
else
	intens = load(inten,'-mat');
end
if isstruct(roi)
	rois = roi; clear roi;
else
	rois = load(roi,'-mat');
end

if isempty(varargin)
	options = struct;
else
	options = varargin{1};
    if length(varargin) > 1,
        n_repeats = varargin{2};
    end
end

options = default(options,'control',[]);
options = default(options,'roitags', []);
options = default(options,'step', true);
options = default(options, 'filter_intens', true);
options = default(options, 'color_order', true);

%% TODO: make these members of an options structure (and give reasonable defaults)

plotdims = [3 3]; % TODO this needs to have a more flexible format SOON
pcutoff = 0.01;
powcutoff = 0.1 ;
minstd = 0.001;
stimnames = rois.header.stim_labels;
nstacks = sum(rois.header.nstacks);
stimons = find(diff(rois.header.stim_lookup(1:sum(rois.header.nstacks)))>0)+1;
stimids = rois.header.stim_lookup(stimons);
% stimons(1) = []; stimids(1) = [];  % to account for first-trial rundown
first_stimuli = length(rois.header.nframes);
edges = [0 (rois.header.nframes./rois.header.frames_per_stack)];
edges = edges(1:end-1);
for i = 1:first_stimuli;
    stimoff(i) = find((stimons - edges(i))>1,1);
end
stimoff = stimons(stimoff);
plotorder = 1:length(stimnames);
% if applying filter
if options.filter_intens
    [b,a] = cheby2(4,25,[1/1000 1/10.1]/(0.5/5));
    tmp1 = [intens.intensities(:,end:-1:1) intens.intensities];
    tmp = filter(b,a,tmp1,[],2);
    intens.intensities = tmp(:,size(intens.intensities,2)+1:end);
end
%% 
mfig = figure('color',[0 0 0], 'position',[10 10 600 800],'name','VOL_ROI_STEPPER: PRESS Q TO QUIT','inverthardcopy','off','paperpositionmode','auto');

if isempty(options.roitags)
	z = 1;
elseif isnumeric(options.roitags{1})
	z = find(options.roitags{1}==[rois.roi_defs.label]);
elseif iscellstr(options.roitags)
	z = strmatch(options.roitags{1},{rois.roi_defs.label},'exact');
end

if isempty(options.roitags)
	roiIdx = [];
else
	roiIdx = getIdx(rois,options.roitags);
end  

if length(roiIdx) > 1
	options.step = false; % TODO: can this be allowed within the stepper loop?
end

if isnumeric(intens.tags{1})
	usz = unique([intens.tags{:}]);
% 	usz = num2str(usz);
else
	usz = unique(intens.tags);
end
sz = length(usz);
xticks = -5:9;
%%
while ~isnan(z)
    delete(get(mfig,'children'));
    thesetraces = nan(length(stimnames), n_repeats, length(xticks));  %optional line, necessary for uneven trials to avoid a final zero line
		for j = 1:min([length(stimnames) prod(plotdims)])
			jor = plotorder(j);
			ax(jor) = subplot(plotdims(1), plotdims(2), j);
			theseon = stimons(stimids==jor);
            theseon = theseon(theseon<(nstacks-9));
			for k = 1:length(theseon)
                if ismember(theseon(k), stimoff)
                    thesetraces(j,k,:) = NaN;
                elseif ~isempty(roiIdx)
					ratio = [rois.roi_defs(roiIdx).roiVolumeInUm3]/sum([rois.roi_defs(roiIdx).roiVolumeInUm3]);
					tmp = intens.intensities(roiIdx,theseon(k)-5:theseon(k)+9);
					tmp = bsxfun(@times,tmp,ratio');
					thesetraces(j,k,:) = mean(tmp,1);  %Here we want NaN for mean comprising weighted points when there is a NaN element 
				else
					if isnumeric(usz)
						thesetraces(j,k,:) = mean(intens.intensities(usz(z)==[intens.tags{:}],theseon(k)-5:theseon(k)+9),1);
					else
						thesetraces(j,k,:) =  mean(intens.intensities(strmatch(usz{z},intens.tags),theseon(k)-5:theseon(k)+9),1);
					end
				end
				
				thisbase(j,k) = nanmean(thesetraces(j,k,1:5));
				thesetraces(j,k,:) = (thesetraces(j,k,:)-thisbase(j,k))/thisbase(j,k);
				thisnormbase(j,k) = nanmean(thesetraces(j,k,1:5));
			end
			% cut out first Ringer's trial (where rundown is there)
			if ~isempty(strmatch(stimnames{jor},'Ringer''s'))
				thesetraces(j,1,:) = NaN;
				thisbase(j,1) = NaN;
				thisnormbase(j,1) = NaN;
			end
			%
			thismean(j,:) = nanmean(squeeze(thesetraces(j,:,:)),1);
			meanpow(j) = nanmean(thismean(j,6:10));
			thisz(j,:) = thismean(j,:)./(minstd+nanstd(squeeze(thesetraces(j,:,:)),[],1));
			thisz(j,isnan(thisz(j,:))) = 0;
			maxz(j) = nanmax(thisz(j,6:10));
			tmp = sort(thisz(j,6:10),'descend');
			meanz(j) = nanmean(tmp(1:end));
            if options.color_order
               set(gcf, 'DefaultAxesColorOrder', cool(size(thesetraces,2)));
%                pre = zeros(3);pre(:,1) = 1;post = zeros(6,3);post(:,3) = 1; conditions = [pre; post]; set(gcf, 'DefaultAxesColorOrder', conditions);
                plot(squeeze(100*thesetraces(j,:,:))','linewidth',2);
            else
                plot(squeeze(100*thesetraces(j,:,:))','linewidth',2,'color',[0.6 0.6 0.6]);
            end
			hold on;
            plot(100*thismean(j,:),'w-','linewidth',2, 'color', [1 1 1]); 
			%plot(100*thismean(j,:),'w-','linewidth',2.5, 'color', [0 0 0]);  %switch color back
			panmax(j) = 100*max(max(thesetraces(j,:,:),[],3),[],2);
			panmin(j) = 100*min(min(thesetraces(j,:,:),[],3),[],2);
			title(stimnames{jor});
			set(get(gca,'title'),'color',[0 0 0]);
			set(ax(jor),'color',[0 0 0],'xcolor',[1 1 1],'ycolor',[1 1 1],'box','off','xtick', 1:5:16, 'xticklabel', {'-5','0','5','10'});
            %set(ax(jor),'color',[1 1 1],'xcolor',[0 0 0],'ycolor',[0 0 0],'box','off','xtick', 1:5:16, 'xticklabel', {'-5','0','5','10'});
		end
		if max(panmax) > 0.1 || max(panmax) < 0.03
			set(ax,'ylim',[min(panmin) max(panmax)]);
    else
        set(ax,'ylim',[-0.08 0.08]);
    end
    % apply a statistical filter to find significant responses
    if ~isempty(options.control)
        if ischar(options.control)
            ctrlidx = strmatch(options.control,stimnames);
        elseif isnumeric(options.control)
            ctrlidx = options.control;
        end
        ctrlidx = find(plotorder==ctrlidx);
        ctrlsamps = squeeze(thesetraces(ctrlidx,2:end,6:10));
        ctrlsz = size(ctrlsamps,1);
        ctrlgroup = repmat({'Control'},1,ctrlsz);
        for i = 1:length(stimnames)
            if i == ctrlidx; testsamps{i} = []; p(i) = 1; continue; end;
            testsamps{i} = squeeze(thesetraces(i,:,6:10));
            testsz{i} = size(testsamps{i},1);
            testgroup{i} = repmat(stimnames(plotorder(i)),1,testsz{i});
            q = cat(1,ctrlsamps,testsamps{i});
            g = cat(2,ctrlgroup,testgroup{i});
            p(i) = kruskalwallis(q',g,'off');
            % if this is significant, make it a color!
            if p(i) <= pcutoff && meanpow(i) >= powcutoff
                tmpline = get(ax(plotorder(i)),'children');
                set(tmpline(2:end),'color',[0.6 0.3 0.3]);
                set(tmpline(1),'color',[1 0 0]);
            end
        end
        
        htax = axes('position',[.4 .95 .2 .05]);
        set(htax,'xtick',[],'ytick',[],'color',[0 0 0],'box','off');
        htitle = text('string',['ROI # ' num2str(z)],'position', [.4 .4],'visible','on','color',[1 1 1]);
    end
    % panel with heatmap
    hmap = axes('position', [0.92 0.3 0.05 .4],'box','on','xtick',[],'ytick',[]);
    %imagesc(maxz([1:12 14 15 13])');
    imagesc(meanz(1:length(stimnames))');
    set(gca,'clim',[-3 3]);
    colormap(colorize_asymmetric([-3 3]));
    
    % panel with legend
    lmap = axes('position',[0.4 0.01 0.2 .05],'box','off','xtick',[],'ytick',[],'color',[0 0 0]);
    set(lmap,'clim',get(hmap,'clim'));
    colorbar('peer',lmap,'north');
    
		if ~options.step
			z = NaN;
		else
			z = keystepper(1:sz,z,struct('q_only',1));
		end
end

function idx = getIdx(roi,lbl)
tags = {roi.roi_defs.label};
if ischar(lbl)
	idx = strmatch(lbl,tags,'exact');
elseif iscellstr(lbl)
	for i = 1:length(lbl)
		idx(i) = strmatch(lbl{i},tags,'exact');
	end
elseif iscell(lbl)
	for i = 1:length(lbl)
		idx(i) = find(lbl{i}==[roi.roi_defs.label]);
	end
end