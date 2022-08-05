function pyramids_imagine(smm, params)
% PYRAMIDS_IMAGINE: create smoothed, subsampled versions of IMAGINE stacks
% Syntax:
%   pyramids_imagine(smm, params)
% where
%   smm is the memory-mapped imagine file (see STACKMM);
%   params may contain the following fields:
%     stacknumbers: a vector of stack numbers to be processed (required);
%     outputpath (default: current directory): a string
%     decimation_schedule: a matrix specifying the degree of
%       smoothing/decimation to be applied, recursively, at
%       each level of the pyramid. The default is
%           4 4 1
%           2 2 1
%           2 2 2
%       which results in a pyramid with 3 levels, the finest level
%       binned in 4-by-4-by-1 (x-y-z) bins, the next an additional
%       2-by-2-by-1 binning (for a total of 8-by-8-by-1), and the
%       coarsest level an additional 2-by-2-by-2 (for a total of
%       16-by-16-by-2).
% Upon output, the directory specified by outputpath will have new
% directories named after the imagine file, with _4x4x1 (etc.) appended.
% Each directory will be populated by .mat files with filenames like
% "stack002.mat", each containing a single image.
%
% See also: IMREDUCE, STACKMM.
  
% Copyright 2006 by Timothy E. Holy
  
  [imaginefile_path,imaginefile_base] = fileparts(smm.filename);
  if isempty(imaginefile_path)
      imaginefile_path = pwd;
  end
  if ~isfield(params,'outputpath') || isempty(params.outputpath)
    params.outputpath = pwd;
  end
  if ~isfield(params,'decimation_schedule')
    params.decimation_schedule = [4 4 1; 2 2 1; 2 2 2];
  end
  if ~exist(params.outputpath,'dir')
    [status,message] = mkdir(params.outputpath);
    if (status)
      warning([params.outputpath ' created.'])
    else
      error(message)
    end
  end
  total_decimation_schedule = cumprod(params.decimation_schedule,1);
  nlevels = size(params.decimation_schedule,1);
  dirnames = cell(1,nlevels);
  for i = 1:nlevels
    dirnames{i} = [params.outputpath filesep imaginefile_base '_pyramids' ...
      num2str(total_decimation_schedule(i,1)) 'x'...
      num2str(total_decimation_schedule(i,2)) 'x'...
      num2str(total_decimation_schedule(i,3))];
    if ~exist(dirnames{i},'dir')
      [status,message] = mkdir(dirnames{i});
      if (status)
        warning([dirnames{i} ' created.'])
      else
        error(message)
      end
    end
  end
  nstacks = length(params.stacknumbers);
  ndigits = ceil(log10(max(params.stacknumbers)));
  filenames = cell(1,nstacks);
  for i = 1:nstacks
    filenames{i} = sprintf(['stack%0' num2str(ndigits) 'd'], ...
			   params.stacknumbers(i));
  end
  h = smm.header;
  spacing0 = [h.um_per_pixel_xy([1 1]) diff(h.piezo_start_stop)/h.frames_per_stack];
  
  for i = 1:nstacks
    im = single(smm(:,:,:,params.stacknumbers(i)));
    for j = 1:nlevels
      im = imreduce(im,params.decimation_schedule(j,:));
      spacing = spacing0 .* total_decimation_schedule(j,:);
      save([dirnames{j} filesep filenames{i}],'im','spacing');
    end
  end
  
