function [roi_defs, ops, N, roi_img] = roi_by_imflow(varargin)
% creates regions of interest in a .imagine movie based on activity
%
% Syntax: [roi_defs, ops, N, roi_img] = roi_by_imflow(smm_in, options)
%                                ...  = roi_by_imflow('previous_save.mat')
%
% Inputs: smm_in: a stackmm object or a string identifying a '.imagine' file
%         options: a structure with the following subfields:
%                 .func: a function handle, either 'imflow_dotproduct_mex'or 'imflow'
%                 .metric: a string specifying either 'dfof' or 'zscore' (default 'dfof')
%                 .dfof_thresh: a scalar specifying the activity dfof threshold (default [] )
%                          (if BOTH dfof_thresh & zcore_thresh left blank, a figure will pop up to help you choose the threshold)
%                 .zscore_thresh: a scalar specifying the activity z threshold (default [] )
%                          (if BOTH dfof_thresh & zcore_thresh are left blank, a figure will pop up to help you choose the threshold)
%                 .sigma: value of the gaussian filter to be applied to the dataset BEFORE thresholding (default [] )
%                 .max_size: an integer value of the maximum dimensions of the evaluated stack (default [] )
%                 .bg_stks: scalar vector specifying stacks to use for background (referenced to stim. onset)
%                              (default -5:-1)
%                 .stim_stks: scalar vecor specifying stacks to use post-stimulus (ref'd to stim. onset)
%                              (default 0:4)
%                 .stim_to_use: scalar vector or cell array of strings specifying stimuli to use
%                              (default 1:length(smm_header.stim_labels) )
%                              (e.g. 1:12, [1 3 5], or {'male urine', 'female urine'} )
%                 .trials_to_use: scalar vector of repeat #s to use (default 1:n_repeats)
%                 .output_file: a filename to write the .roidef info to (default [smm_header.filename '.imflow.roidef']);
%                 .std_offset: a scalar indicating the offset to be added to the denominator when calculating zscore (i.e. mean(dfof)/(std(dfof)+std_offset))
%                               (default 1e-3) 
%                 .bad_stacks: an integer vector containing valid stack indices which should be ignored
%                 .isbad:  logical array from structure output find_bad_frames.m, converted to .bad_stacks; 
%                               note that after registration these values 
%                               have usually become part of the new .cam
%                               movie and do not necessarily need to be
%                               exluded
%                 .min_pixels: integer minimum ROI pixel size.  Automatic ROIs with fewer than this number will be ignored (deleted). default: 10
%                 .edge_mask: a vector containing integer values to mask on
%                   each edge [bottom top left right].
%                 
% Outputs: roi_defs: a structure variable with the following subfields:
%                 .label: the integer ID of the ROI
%                 .centerInPixels: the pixel-wise center of mass of the ROI
%                 .vtxInPixels: the pixel-wise vertices of the shell around the ROI
%                 .centerInUm: the micron-valued center of mass of the ROI
%                 .vtxInUm: the micron-valued center of the shell around the ROI
%                 .roiVolumeInUm3: the volume of the ROI in cubic microns
%                 .pixels: just the indices of the pixels in the ROI (@ full frame)
%                 .weight: a matrix of same size as .pixels specifying the "weight" to be given to each pixel
%                          when computing ROI sum activity
%          ops: the options structure ultimately used to create the ROIs (i.e. with defaults/choices filled in)
%          N: the number of ROIs identified (i.e. length(roi_defs))
%          roi_img: an image of the ROIs, using the same dimensions as the 
%                   image evaluated for ROIs (may be lower-res than .imagine movie)
%
% See also imflow, imflow_dotproduct, stackmm

%
% Copyright 2011 Julian P Meeks (Timothy Holy Laboratory)
% History:
%   2011_08_07: wrote it (JPM)
%
%

% check args for options
if nargin > 2
	warning('more than two arguments supplied, only the first 2 are used');
elseif nargin == 2
	if ~isstruct(varargin{2})
		error('second argument must be a structure variable, see help for details');
	else
		options = varargin{2};
	end
else
	options = struct; 
end

options = default(options, 'func', 'imflow_dotproduct');
options = default(options, 'metric', 'dfof');
options = default(options, 'dfof_thresh', [] );
options = default(options, 'zscore_thresh', [] );
options = default(options, 'sigma', [] );

options = default(options, 'bg_stks', -5:-1);
options = default(options, 'stim_stks', 0:4);
options = default(options, 'stim_to_use', [] );
options = default(options, 'trials_to_use',[] );
options = default(options, 'output_file',[]);
options = default(options, 'std_offset',1e-3);
options = default(options, 'bad_stacks',[]);
options = default(options, 'isbad',[]);
options = default(options, 'min_pixels', 10);
options = default(options, 'edge_mask', [0 0 0 0]);

% check stack input
if nargin < 1
	smm_string = uigetfile('*.imagine', 'Please select an input file');
	smm_in = stackmm(smm_string);
	smm_header = smm_in.header;
	contin = false;
elseif isstr(varargin{1})
	if isempty(strfind(varargin{1},'.imagine'))
		new_ops = options;
		load(varargin{1});
		if nargin > 1
			warning('Please note: when specifying a.mat file, options structure is overwritten');
		end
%     fldnames = fieldnames(options);
% 		for i = 1:fieldnames
% 			if isfield(new_ops,fldnames{i})
% 				options.(fldnames{i}) = new_ops.(fldnames{i});
% 			end
% 		end
		contin = true;
	else
		smm_string = varargin{1};
		smm_in = stackmm(smm_string);
		smm_header = smm_in.header;
		contin = false;
	end
elseif isobject(varargin{1})
	smm_in = varargin{1};
  smm_header = smm_in.header;  %%gfh
	contin = false;
else
	error('unable to parse stack input. See help.');
end

% extract important info from smm_in.header
if ~contin  % THIS OPTION ALLOWS A USER TO SAVE and RETURN TO PREV. SAVED DATA
	
smm_sz = smm_in.size;
n_stacks = smm_sz(4);
im_dims = smm_sz(1:3);
n_stim = length(smm_header.stim_labels);
s_lab = smm_header.stim_labels;
if isempty(options.stim_to_use)
	options.stim_to_use = 1:length(s_lab);
	stu = options.stim_to_use;
elseif iscellstr(options.stim_to_use)
	for i = 1:n_stim
		stu(i) = strmatch(options.stim_to_use{i}, s_lab,'exact');
	end
else
	stu = options.stim_to_use;
end

options = default(options, 'max_size', im_dims);

nstims_chosen = length(stu);
stims = smm_header.stim_lookup;
stims = stims(1:smm_sz(4));
onsets = find(diff(stims)>0)+1;
stim_order = stims(onsets);
if isfield(smm_header, 'trial_lookup')
  trialnum = smm_header.trial_lookup;
  trialnum = trialnum(onsets);
  n_repeats = max(trialnum);
else
  [sorted_stim_order, sorti] = sort(stim_order);
  edges = find(diff(sorted_stim_order)==1)'; 
  edges = [0 edges size(stim_order,1)];
  stimulus = struct('trials',{});
  for i = 1:size(stu,2)
    by_trial = cell(1,(edges(i+1)-edges(i)));
    by_trial = onsets(sorti((edges(i)+1):edges(i+1)));
    stimulus(stu(i)).trials = deal(by_trial');
  end
  n_repeats = max(arrayfun(@(x)(length(x.trials)),stimulus));
end

    
if ~isempty(options.trials_to_use)
	ttu = options.trials_to_use;
else
	ttu = 1:n_repeats;
end
ntrials_chosen = length(ttu);

% choose stack indices for background, stimulated
bkgnd = zeros(nstims_chosen,ntrials_chosen,length(options.bg_stks));
stimd = zeros(nstims_chosen,ntrials_chosen,length(options.stim_stks));
for i = 1:nstims_chosen
	stmp = stim_order==stu(i);
	otmp = onsets(stmp);
	for j = 1:ntrials_chosen
		bkgnd(i,j,:) = otmp(ttu(j))+options.bg_stks;
		stimd(i,j,:) = otmp(ttu(j))+options.stim_stks;
	end
end; clear stmp otmp;

% convert isbad to bad_stacks
if ~isempty(options.isbad)
  options.bad_stacks = unique(options.bad_stacks, find(sum(options.isbad,1)>0));
end

% set decimation schedule, if any
stk = [];   % placeholder
if ~isempty(options.max_size)
	maxsz = options.max_size;
else	
	maxsz = im_dims;
end
restrict_sched = [0 0 0];
tmp = im_dims;
while any(tmp > maxsz)
	restrict_sched = restrict_sched + (tmp > maxsz);
	tmp(tmp>maxsz) = ceil(tmp(tmp>maxsz)/2);
end; clear tmp;
imdec = @(stk)do_restrict(stk,restrict_sched);

% set filter plans, if any
if ~isempty(options.sigma)
	filt = options.sigma;
else
	filt = [0 0 0];
end
imfilt = @(stk)do_filter(stk,filt);

% calculate image array
im2flow.dfof = single(zeros([size(imdec(single(smm_in(:,:,:,1)))) nstims_chosen]));
im2flow.zscore = im2flow.dfof;
bg = single(zeros(im_dims));
act = bg;
dfof = single(zeros([size(im2flow.dfof(:,:,:,1)) ntrials_chosen]));
tic;
prog  = progress(struct('progress',0,'max',nstims_chosen*ntrials_chosen));
for i = 1:nstims_chosen
	for j = 1:ntrials_chosen
		thesebk = squeeze(bkgnd(i,j,:)); thesebk([find(thesebk==1) find(thesebk==2)])=[]; % note: first stacks are almost always oversaturated or bad
    [~, to_blank] = intersect(thesebk,options.bad_stacks); 
    thesebk(to_blank)=[]; 
    if isempty(thesebk)
      warning(['no background stacks to average for stim ' num2str(i) ' trial ' num2str(j) ' starting at frame ' num2str(bkgnd(i,j,1)) '.\n']);
			dfof(:,:,:,j) = NaN;
			continue;
    end
		thesestm = squeeze(stimd(i,j,:));
    [~, to_blank] = intersect(thesestm,options.bad_stacks); 
    thesestm(to_blank)=[]; 
    if isempty(thesestm)
      warning(['no stim stacks to average for stim ' num2str(i) ' trial ' num2str(j) ' starting at frame ' num2str(stimd(i,j,1)) '.\n']);
			dfof(:,:,:,j) = NaN;
			continue;
		end
		bg = nanmean(single(smm_in(:,:,:,thesebk)),4); zro = bg==0;%bg(bg==0)=NaN;
	  bias(i,j) = min(min(min(nonzeros(bg),[],3),[],2),[],1);  % bias indexed to keep track during testing
		bg = bg - bias(i,j); bg(bg<0) = 0;
		act = nanmean(single(smm_in(:,:,:,thesestm)),4); %act(act==0)=NaN;
		act = act-bias(i,j); act(act<0)=0;
		tmpdfof = (act-bg)./bg; % this is dfof
		tmpdfof(isinf(tmpdfof))=0;
		zro = floor(imdec(single(zro))); zro(isnan(zro))=0; zro = logical(zro);
		tmpdfof = imfilt(imdec(tmpdfof));
		tmpdfof(zro) = NaN;
		dfof(:,:,:,j) = tmpdfof;
		prog.progress = prog.progress+1;
		prog = progress(prog);
	end
	mean_dfof = nanmean(dfof,4);
	stdev_dfof = nanstd(dfof,[],4);
  zscore_dfof = mean_dfof./(stdev_dfof+options.std_offset); zscore_dfof(isinf(zscore_dfof)) = 0;
	im2flow.zscore(:,:,:,i) = zscore_dfof;
	im2flow.dfof(:,:,:,i) = mean_dfof;		
	if isempty(options.metric)
		options.metric = 'dfof';
	end
end
prog.progress = -1;
prog = progress(prog); if ishandle(prog.handle); delete(prog.handle); end; clear prog;
toc;
clear bg act dfof mean_dfof stdev_dfof zscore_dfof varargin;
save(['roi_by_imflow_temp.' options.metric '.' datestr(now,'yyyymmddHHMM') '.mat'],'-v7.3');

end % END OF (if ~contin) SECTION

% set threshold (support for user if they have supplied no threshold)
if isempty(options.dfof_thresh) || isempty(options.zscore_thresh)
	new_options = threshold_stack_gui(im2flow,options,s_lab);
	if isempty(new_options)
		warning('you did not choose a threshold, exiting program.  if you wish to start the GUI over with the same settings, please supply the ''roi_by_imflow_temp'' filename in your local directory in place of the .imagine filename');
		roi_defs = 0;
		N = 0;
		ops = options;
		roi_img = [];
		return;
	elseif isempty(new_options.dfof_thresh) || isempty(new_options.zscore_thresh)
		warning('you did not choose a threshold, exiting program.  if you wish to start the GUI over with the same settings, please supply the ''roi_by_imflow_temp'' filename in your local directory in place of the .imagine filename');
	  roi_defs = 0;
		N = 0;
		roi_img = [];
		return;
	end
	dfof_thresh = new_options.dfof_thresh;
	zscore_thresh = new_options.zscore_thresh;
    edge_mask = new_options.edge_mask;
	options.dfof_thresh = dfof_thresh;
	options.zscore_thresh = zscore_thresh;
    options.edge_mask = edge_mask;
    
else
	if isempty(options.dfof_thresh)
		dfof_thresh = 0;
		options.dfof_thresh = dfof_thresh;
	else
	  dfof_thresh = options.dfof_thresh;
	end
	if isempty(options.zscore_thresh)
		zscore_thresh = 0;
		options.dfof_thresh = zscore_thresh;
	else
	  zscore_thresh = options.zscore_thresh;
    end
    if ~isfield(options, 'edge_mask')  % This is required for reloading old files.
        options.edge_mask = [0 0 0 0];
    end
end

pause(0.1);

if ~isfield(options,'metric')
	if ~isempty(options.metric)
		metric = options.metric;
	else
		metric = 'zscore';
	end
else
	metric = 'zscore';
end
metric = questdlg('Which metric should be used for calculating weights?','Choose weight metric','zscore','dfof',metric);
options.metric = metric; clear metric;

i2f = im2flow.(options.metric);
zro = isnan(max(im2flow.dfof,[],4));
i2f(repmat(zro,[1 1 1 size(i2f,4)])) = 0;
zro = isnan(max(im2flow.zscore,[],4));
i2f(repmat(zro,[1 1 1 size(i2f,4)])) = 0;
maskEdges = false(maxsz(1), maxsz(2));
maskEdges(options.edge_mask(2):(maxsz(1)-options.edge_mask(3)),...
    options.edge_mask(3):(maxsz(2)-options.edge_mask(4)))=true;
	
% apply any thresholds
if length(zscore_thresh) == 1 && length(dfof_thresh) == 1
	mi2f = max(im2flow.dfof,[],4);
	subthresh = mi2f < dfof_thresh;
	i2f(repmat(subthresh,[1 1 1 size(i2f,4)])) = 0; clear mi2f subthresh;
	mi2f = max(im2flow.zscore,[],4);
	subthresh = mi2f < zscore_thresh;
	i2f(repmat(subthresh,[1 1 1 size(i2f,4)])) = 0; clear mi2f subthresh;
    i2f = bsxfun(@times,i2f,maskEdges);
elseif size(zscore_thresh,1)==smm_sz(3) && size(zscore_thresh,1)==smm_sz(3)
	for z = 1:smm_sz(3)
		multi = [];
		for i = 1:length(options.stim_to_use)
			dmulti = im2flow.dfof(:,:,z,i) > dfof_thresh(z,i);
			zmulti = im2flow.zscore(:,:,z,i) > zscore_thresh(z,i); % end after this line?
			if isempty(multi)
				multi = dmulti.*zmulti;
			else
				multi = multi | dmulti.*zmulti;                
            end
            multi = multi.*maskEdges; %probably error, may need bsxfun
		end
		mask(:,:,z) = multi;
        
	end
	i2f = bsxfun(@times,i2f,mask);
else
	error('did not know how to interpret generated zscore!');
	return;
end

if ~isfield(options,'func')
	if ~isempty(options.func)
		func = options.func;
	else
		func = 'imflow_dotproduct';
	end
else
	func = 'imflow_dotproduct';
end
func = questdlg('Which function should be used for flowing pixels?','Choose flow function','imflow_dotproduct','imflow',func);
options.func = func; clear func;

% use the imflow function
switch options.func
	case 'imflow_dotproduct'
		temp = permute(i2f,[4 1 2 3]);
		map = single(imflow_dotproduct_mex(temp));
		map = killedges(map,0);
		zro = map==0;
		map(zro)=find(zro);
	case 'imflow'
		if length(size(i2f)) > 3
			temp = single(squeeze(max(i2f,[],4)));
			%temp(temp==0) = NaN;
			map = single(imflow_mex(temp));
			%map = map+1+prod(size(map(:,:,1)));
 			map = killedges(map,0);
			zro = map==0;
			map(zro)=find(zro);
		end
end;

% flow the pixels to their max
mapOld = []; 
while ~isequal(mapOld,map)
	mapOld = map;
	map = map(map);
end
clear temp;
I = squeeze(max(i2f,[],4));
I(zro)=NaN;
I = I(map);  roi_img = I; % OUTPUT 4

% identify number ROIs
this = 0;
count = 0;
RoI = cell(1);
useful = find(I>0);
this = find(map(useful)==useful);
pks = useful(this);
N = length(pks);

for i = 1:N
	n_to_pk(i) = length(find(map(useful)==pks(i)));
end

% exclude small glomeruli as per options.min_pixels
options = default(options,'min_pixels',10); % since early users did not have this option;
this(n_to_pk < options.min_pixels) = [];
pks(n_to_pk < options.min_pixels) = [];
n_to_pk(n_to_pk < options.min_pixels) = [];
N = length(pks); % OUTPUT 3

prog  = progress(struct('progress',0,'max',N+10));

[~, pksort] = sort(n_to_pk,'descend');

for i = 1:N
	thispk = pks(pksort(i));
	RoI{i} = find(map==thispk);
end

% put together roi_defs structure
roi_defs = struct;
eek = array_prolong(i2f,[im_dims size(i2f,4)]);

prog.progress = 10;
prog = progress(prog);

for i = 1:N
	vtx = [];
	dummy = logical(zeros(size(I)));
	decimations = ceil(smm_sz(1:3)./size(I));
	pixel_spacing = smm_header.pixel_spacing;
	roi_defs(i).label = num2str(i);
  dummy(RoI{i}') = 1;	
	deci = log2(decimations);
	if any(decimations>0)
		tmp = array_prolong(double(dummy),im_dims);
		tmp = logical(round(tmp));
	end
	pixi = find(tmp);
	[roi_defs(i).pixels(:,1) roi_defs(i).pixels(:,2) roi_defs(i).pixels(:,3)] = ind2sub(im_dims,pixi);
	% since we have the data now, compute the normalized correlation with the mean (for purposes of "weighting" the ROI later)	
	snip = [];
	for index = 1:size(i2f,4)
		cut = eek(:,:,:,index);
		snip(:,index) = cut(pixi);
	end
	rmean = mean(snip,1);
	c = sum(bsxfun(@times,snip,rmean),2);
	c(c<0)=0;
	roi_defs(i).weight = c/sum(c);
	% junk = eek(:,:,:,1); junk = 0; junk(pixi) = c;
	roi_defs(i).centerInPixels = round(mean(roi_defs(i).pixels,1));
	roi_defs(i).centerInUm = roi_defs(i).centerInPixels.*pixel_spacing;
	perim = bwperim(tmp,4);
	% a = single(dummy); a(:) = 0; a(dummy) = -1; a(perim)=1; figure; z = 1; while ~isnan(z); set(gca,'clim',[-1 1]);imagesc(a(:,:,z));colormap(colorize_asymmetric([-1 1])); z = keystepper(1:51,z);end;
	[vtx(:,1), vtx(:,2) vtx(:,3)] = ind2sub(im_dims,find(perim));
	%vtx(:,1) = round(vtx(:,1)*decimations(1));
	%vtx(:,2) = round(vtx(:,2)*decimations(2));
	%vtx(:,3) = round(vtx(:,3)*decimations(3));
	% need to prep vertices to be read by existing infrastructure if possible
	[s, si] = sort(vtx(:,3),'ascend');
	vtx = vtx(si,:);
	for j = 1:im_dims(3)
		roi_defs(i).vtxInPixels{j} = [];
		these = find(s==j);
		if ~isempty(these)
			thesevtx = vtx(these,:);
			while ~isempty(thesevtx)
				this = thesevtx(1,:);
				roi_defs(i).vtxInPixels{j} = [roi_defs(i).vtxInPixels{j}; this];
				thesevtx(1,:) = [];
				dist = sqrdist(this', thesevtx');
				[~, sd] = sort(dist,'ascend');
				thesevtx = thesevtx(sd,:);
			end
		end
		if ~isempty(roi_defs(i).vtxInPixels{j})
			roi_defs(i).vtxInUm{j}(:,1) = roi_defs(i).vtxInPixels{j}(:,1)*pixel_spacing(1);
			roi_defs(i).vtxInUm{j}(:,2) = roi_defs(i).vtxInPixels{j}(:,2)*pixel_spacing(2);
			roi_defs(i).vtxInUm{j}(:,3) = roi_defs(i).vtxInPixels{j}(:,3)*pixel_spacing(3);
		else
			roi_defs(i).vtxInUm{j} = [];
		end
	end	
	roi_defs(i).roiVolumeInUm3 = size(roi_defs(i).pixels,1)*prod(pixel_spacing);
	prog.progress = i+10;
	prog = progress(prog);
end

for i = 1:length(roi_defs)
	roi_defs(i).label = num2str(i);
end

prog.progress = -1;
prog = progress(prog); if ishandle(prog.handle); delete(prog.handle); end; clear prog;

% save .roidef file
if ~isfield(options,'output_file')
	options.output_file = [];
end
if isempty(options.output_file)
	options.output_file = [smm_header.header_filename '.imflow.roidef'];
end
if iscellstr(options.output_file)
	options.output_file = options.output_file{1};
end
if exist(options.output_file,'file')
	response = questdlg('File exists, what do we do? :', 'Existing File Found', 'Overwrite','Choose other name','Choose other name');
	if ~isempty(strmatch('Overwrite',response)) 
		delete(options.output_file);
	elseif ~isempty(strmatch('Choose other name',response))
		[filename path] = uiputfile('*.roidef', 'Choose an output file name');
		options.output_file = [path filename];
	else
		warning('roi_defs not saved to file.  Did you mean to do this?');
		options.output_file = [];
	end
end
header = smm_header;
pixelPerUm = 1./pixel_spacing;
if ~isempty(options.output_file)
	save(options.output_file,'roi_defs','header','pixelPerUm','dfof_thresh','zscore_thresh');
end
	
ops = options; % OUTPUT 2
end

function imout = do_restrict(imin,restrict)
  imout = imin;
	cp = single(imin~=0);
  while any(restrict)
		to_r = restrict>0;
    imout = array_restrict(imout,to_r);
		cp = array_restrict(cp,to_r);
    restrict = restrict-to_r;
	end
	cp(cp>0) = 1;
	imout = imout.*cp;
	imout(imout==0)=NaN;
end

function imout = do_filter(imin,filt)
  n = isnan(imin);
  imout = imfilter_gaussian(imin,filt);
	imout(n)=NaN;
end