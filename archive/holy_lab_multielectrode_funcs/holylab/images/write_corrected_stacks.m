function write_corrected_stacks(smm,filenameout,options)
% write_corrected_stacks: utility for balancing intensity, removing stripes, clims, and cropping
%
% Syntax:
%   write_corrected_stacks(smm,filenameout,options)
% where
%   smm is a stackmm object
%   filenameout is the name of the file you want to write (if there are '.'
%     in the filename, be sure to supply an extension or everything after the
%     '.' will be stripped away)
%   options may have the following fields:
%     stacknums (default 1): the list of stack numbers you wish to process
%       with this function
%     ROI (default whole stack): a 3-element cell array giving the x,y,z
%       coordinates of the region you want to keep, e.g.,
%             options.ROI = {':','20:2065',':'};
%     stripe_mask (default none): if you want to correct stripes, supply a
%       stripe mask. The mask must be created to work on the original
%       frame size (i.e., before applying ROI).
%     custom_function: if supplied, the pixel intensity is converted in
%       the following manner:
%           imout = options.custom_function(im);
%       For example, one could set
%           options.custom_function = log(1+im);
%       to do a log-transform.
%     clim: a 2-element vector specifying the range of intensities. If
%       supplied, the image intensity is truncated to within the clim range.
%       clim(1) is subtracted so that the lowest intensity in the image is
%       zero. clim is applied after all of the above transformations.
%     precision: a string, if supplied the output pixels are converted to
%       this data format (e.g., precision = 'uint16'). If clim is also
%       provided, and if the upper value is larger than this type can
%       support, pixels values are multiplied by the ratio
%       (max_value/diff(clim)) before writing.
%     center_factor (default 0.45): the relative position of the peak
%       intensity of the light sheet
%     width_factor (default 0.28): the relative width of the light sheet

% Copyright 2011 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  options = default(options,...
    'stacknums',1,...
    'ROI',{':',':',':'},...
    'max_frames_simultaneously',100,...
    'center_factor',0.45,...
    'width_factor',0.28,...
    'stripe_mask',[]);
  
  % Determine the dimensions of the cropped region
  sz = smm.size;
  for dimIndex = 1:3
    if ischar(options.ROI{dimIndex})
      options.ROI{dimIndex} = 1:sz(dimIndex);
    end
  end
  imROI = options.ROI(1:2);
  
  % Split the frame index up into chunks of manageable size
  splits = 0:options.max_frames_simultaneously:length(options.ROI{3});
  if (splits(end) < options.ROI{3}(end))
    splits(end+1) = options.ROI{3}(end);
  end
  frameList = cell(1,length(splits)-1);
  for i = 1:length(frameList)
    frameList{i} = options.ROI{3}(splits(i)+1:splits(i+1));
  end
  
  % Generate the lightsheet intensity correction factor
  x = 1:sz(1); x = x(:);
  xc = options.center_factor*sz(1);
  sigma = options.width_factor*sz(1);
  f = sqrt(1+(x-xc).^2/sigma^2);
  
  % Parse the filename
  [pth,basename] = fileparts(filenameout);
  if isempty(pth)
    baseout = basename;
  else
    baseout = [pth filesep basename];
  end
  
  % Open the output raw file
  [fid,msg] = fopen([baseout '.cam'],'w');
  if (fid < 1)
    error(msg)
  end
  closeOnExit = onCleanup(@() fclose(fid));
  
  % Write the .imagine file header
  header = smm.header;
  if isfield(options,'precision')
    header.prec = options.precision;
  end
  header.height = length(options.ROI{2});
  header.width = length(options.ROI{1});
  header.nframes = length(options.ROI{3});
  header = rmfield(header,'depth');
  header.nstacks = length(options.stacknums);
  save(baseout,'header')
  
  % Calculate the scale factor
  defscale = 1;
  if isfield(options,'clim')
    if (header.prec(1) == 'u' || header.prec(1) == 'i')
      mx = intmax(header.prec);
      defscale = double(mx)/diff(options.clim);
    end
  end
  options = default(options,'scalefactor',defscale);
  
  for stackCounter = 1:length(options.stacknums)
    stackIndex = options.stacknums(stackCounter);
    for chunkIndex = 1:length(frameList)
      fprintf('%d%%...',round(100*chunkIndex/length(frameList)));
      im = smm(:,:,frameList{chunkIndex},stackIndex);
      % Correct the intensity
      im = bsxfun(@times,double(im),f);
      % Correct the stripes
      if ~isempty(options.stripe_mask)
        for frameIndex = 1:size(im,3)
          im(:,:,frameIndex) = imfilter_fourier_mask_apply(im(:,:,frameIndex)+1,options.stripe_mask,struct('log',true));
        end
      end
      % Crop the image
      im = im(imROI{:},':');
      % Apply the scaling function
      if isfield(options,'custom_function')
        im = options.custom_function(im);
      end
      % Apply clims
      if isfield(options,'clim')
        im(im > options.clim(2)) = options.clim(2);
        im = im - options.clim(1);
        im(im < 0) = 0;  % for non-uint data types
        % Scale if needed before casting
        if (options.scalefactor ~= 1)
          im = cast(options.scalefactor*im,header.prec);
        end
      end
      % Save to disk
      count = fwrite(fid,im,header.prec);
      if (count < numel(im))
        error('Error writing file')
      end
    end
  end
  fprintf('\nDone.\n');
end

      
        