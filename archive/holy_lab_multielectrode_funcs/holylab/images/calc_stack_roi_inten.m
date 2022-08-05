function [intensities tags] = calc_stack_roi_inten(basename, roiDefFile, options)
%calc_stack_roi_inten takes roi definitions within stack space, created by stack_roi,
%and calculates the mean pixel intensity within each of those regions over time. The
%output is an RxS matrix, where R is the # of rois and S is the # of stacks
%in the experiment.
% Syntax:
%   intensities = calc_stack_roi_inten(basename, roiDefFile)
%   intensities = calc_stack_roi_inten(basename, roiDefFile, options)
% where:
%   basename is the original basic filename from the imagine experimental
%             data (also the header name), this information can either be passed as a
%             primary string (e.g. '2007-07-07-7') or as a named variable that contains this string.
%             BTW:  the header file from the original data should automatically be saved
%             as a .mat or .imagine file into the parent directory for this job by
%             the previous processing function (registration, etc.). If not, you should do this.
%             NOTE: if the previous processing function modified the shape
%             of the data (e.g. dropped the first and last stacks) the
%             copied header needs to be edited to reflect this (you may need to check that
%             this is so if the program barfs at you).
%   roiDefFile is the name of the roi definition file created by stack_roi
%   options is a structure indicating some of the extra behavior you would
%             or would not like calc_stack_roi_inten to exhibit.
%        .show_stack_progress if set to 1(true) causes an indicator of
%             progress through stacks to be indicated on the command screen.
%             Default is true.
%        .show_roi_progress if set to 1(true) causes an indicator of
%             progress through rois to be indicated on the command screen.
%             Default is false.
%        .save_intensity_file if set to 1(true) causes a copy of the RxS
%             outcome matrix to be save to disk as a .mat file. Default is true.
%        .savedir if provided allows you to indicate a different save to
%             directory that the parent data directory (default is parent
%             dir).  This function can be used to save the output matrix to another
%             directory when you do not have write permissions in the parent data
%             directory.
%        .gaussian (defaut is true). This is parsed to function
%             'roimeasure'. If true, use gaussian weight (centered at roi, spread as sigmal =r)
%             to caculate the integrated roi intensity from an extended 8r-by-8r square region. 
%             If false, simply caculate the average intensity from the roi region.
%
%Copywrite 2006 by Jason Guo & Terrence Holekamp

% Revision History
% Modification 11/2/2010 (Julian P Meeks)
%  added second output variable (tags) to hold "tag" information
%  corresponding to each ROI in "intensities".  This will be of
%  dimensions 1 x size(intensities,1) (or the # of rois)
%  This is to allow straightforward analysis of ROIs spanning multiple
%  stacks.

if(nargin==2)
    options.show_stack_progress=1;
    options.show_roi_progress=0;
    options.save_intensity_file=1;
    options.savedir = cd;
    options.gaussian = 1;
end
if(~isfield(options, 'show_stack_progress'))
    options.show_stack_progress=1;
end
if(~isfield(options, 'show_roi_progress'))
    options.show_roi_progress=0;
end
if(~isfield(options, 'save_intensity_file'))
    options.save_intensity_file=1;
end
if(~isfield(options, 'savedir'))
    options.savedir = cd;
end
if(~isfield(options,'gaussian'))
    options.gaussian = 1;
end

tt=load(roiDefFile, '-mat');
header=tt.header;
roi_defs=tt.roi_defs;
pixelPerUm=tt.pixelPerUm;
do_tform = 0;
if isfield(tt,'tform_info')
  do_tform = 1;
  tform_info=tt.tform_info;
end
  


% to warn user when headerFile is not same as header.header_filename
%    and re-read header from headerFile
if(~isequal(basename, header.header_filename))
    warning('input header file may not be same as the header used to define ROIs');
    header=imreadheader(basename);
end

% TODO: check if deconvolution is used, if used, recalc z position is pixel

stack = stackmm(basename);
stackSize=stack.size;
if isfield(header, 'compound_cam')
  header.nstacks = sum(header.nstacks);
end
intensities=ones(length(roi_defs), header.nstacks)*inf;

for idxStack=1:header.nstacks
    clear imageDatas; % force matlab to free mem
    imageDatas = cell(1,stackSize(3)); %pre-allocate for speed
    for idxFrame=1:stackSize(3)
        imageDatas{idxFrame}=squeeze(stack(:,:,idxFrame,idxStack));
        imageDatas{idxFrame}=imageDatas{idxFrame}';
    end

    if(options.show_stack_progress)
        % TODO: use fprintf
        disp(['now processing stack ' num2str(idxStack) ' of ' num2str(header.nstacks) ' from ' basename]);
    end
    for idxRoi=1:length(roi_defs)
        if(options.show_roi_progress && mod(idxRoi,10)==1)
            % TODO: use fprintf
            disp(['now processing roi ' num2str(idxRoi) ' of ' num2str(length(roi_defs))]);
        end

        cur_roi_def=roi_defs(idxRoi);
        if do_tform
          tform=interpolate3dTform(tform_info, idxStack-1);
          cur_roi_def.posInPixel=tformfwd(tform, cur_roi_def.posInPixel);
          cur_roi_def.posInPixel(3)=max(min(round(cur_roi_def.posInPixel(3)), stackSize(3)-1), 0);
        end
        
        if isfield(cur_roi_def,'posInPixel');
            intensities(idxRoi, idxStack)=roimeasure(imageDatas{cur_roi_def.posInPixel(3)}, convertRoiFormat(cur_roi_def),[],options);
				elseif isfield(cur_roi_def,'pixels')
					  roiFmt = convertRoiFormat(cur_roi_def);
						roiFmt.size = stackSize(1:3);
						intensities(idxRoi,idxStack)=0;
						intensities(idxRoi,idxStack) = roimeasure(imageDatas,roiFmt);
				elseif isfield(cur_roi_def,'centerInPixels')
            roiFmt = convertRoiFormat(cur_roi_def);
            intensities(idxRoi, idxStack) = 0;
            for idxFrames = roiFmt.keep
                tmpRoi = roiFmt; tmpRoi.vtxInPixels = tmpRoi.vtxInPixels{idxFrames};
                intensities(idxRoi, idxStack)=intensities(idxRoi, idxStack)+...
                    roimeasure(imageDatas{idxFrames}, tmpRoi);
						end
        else
            error('.roidef file is missing pixel position information\n');
        end

    end % for, each roi
end % for, each stack

tags = zeros([1 length(tt.roi_defs)]);
for i = 1:length(tt.roi_defs)
    tags(i) = tt.roi_defs(i).label;
end

if options.save_intensity_file
    c = cd;
    cd (options.savedir)
    save([basename '_roi_intensities.intensity'],'intensities','tags')
    cd(c);
end


function result=convertRoiFormat(roi_def)
if isfield(roi_def,'posInPixel')
    result.type='c';
    result.label=roi_def.label;
    result.x=roi_def.posInPixel(1);
    result.y=roi_def.posInPixel(2);
    result.xyradius=roi_def.xyradiusInPixel;
elseif isfield(roi_def,'pixels')
	  result.type='p';
		result.label = roi_def.label;
		result.pixels = roi_def.pixels;
		result.keep = [];
		for i = 1:length(roi_def.vtxInPixels);
			  if ~isempty(roi_def.vtxInPixels{i})
				  result.keep = [result.keep i];
			  end
		end
		result.weight = roi_def.weight;
elseif isfield(roi_def,'centerInPixels')
    result.type='v';
    result.label=roi_def.label;
    result.vtxInPixels=roi_def.vtxInPixels;
    result.keep = [];
    for i = 1:length(roi_def.vtxInPixels)
        if ~isempty(roi_def.vtxInPixels{i})
            result.keep = [result.keep i];
        end
		end
else
    error('.roidef file is missing pixel position information\n');
end
% TODO: necessary?

result.tform=unit_tform(2);




