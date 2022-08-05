function imagine2med(imaginefile,ranges,basenameout)
% IMAGINE2MED: convert imagine files to medical image file formats
%
% Currently supported:  Mayo Analyze (.img), NIfTI (.nii), and DICOM
% (.dcm) files.
% Limitations: it seems that the Matlab DICOM writer can't handle
% stacks---it separates them out to single image planes, and that seems
% awkward. So, .nii is recommended.  Also, for now the DICOM files don't
% have information about physical coordinates
%
% Syntax:
%   imagine2med(imaginefile,ranges,basenameout)
% where
%   imaginefile is the name of the .imagine file
%   ranges is a cell array containing the selection of voxels to export.
%     For example, if you wanted to export stack numbers 2-10 in their
%     entirety, then
%       ranges = {':',':',':',2:10}
%   basenameout is the base filename, optionally with extension, that you
%   want to save to. In particular,
%      basenameout = ~/imgdata/export.nii
%   will result in stacks written with the names
%      export02.nii, export03.nii, etc.
% 
% This makes use of the SPM toolbox---see
% http://www.fil.ion.ucl.ac.uk/spm/software/spm5/
  
% Copyright 2006 by Timothy E. Holy
  
  % Work around a bug in the SPM code that doesn't allow you to write to
  % a different directory
  [pathname,filename,extname] = fileparts(basenameout);
  d_current = pwd;
  cd(pathname);
  if isempty(extname)
    extname = '.nii';  % by default write NIfTI files
  end
  if isempty(strmatch(extname,{'.img','.nii','.dcm'}))
    error(['Can''t write files with extension ' extname]);
  end
  use_spm = false;
  if ~isempty(strmatch(extname,{'.img','.nii'}))
    use_spm = true;
  end
  if (use_spm && ~exist('spm_write_vol','file'))
    addpath ~/matlab/spm
  end

  % Do some parsing of the header file to extract physical dimensions
  % So far, this is only useful when use_spm is true
  smm = stackmm(imaginefile);
  h = smm.header;
  if strcmp(h.machfmt,'n')
    h.machfmt = 'l';
  end
  v.dt = [spm_type(h.prec) strcmp(h.machfmt,'b')];
  sz = smm.size;
  v.dim = sz(1:3);
  v.mat = eye(4);
  v.mat(1,1) = h.um_per_pixel_xy;
  v.mat(2,2) = h.um_per_pixel_xy;
  v.mat(3,3) = abs(diff(h.piezo_start_stop))/(h.depth-1);
  
  sz = smm.size;
  if strcmp(ranges{4},':')
    ranges{4} = 1:sz(end);
  end
  nstacks = length(ranges{4});
  nchars = ceil(log10(max(ranges{4})));  % # chars needed to encode stack#
  pb_options = struct('progress',0,'max',nstacks,'what','Converting stacks...');
  for i = 1:nstacks
    pb_options = progress_bar(pb_options);
    rangestmp = ranges;
    rangestmp{4} = ranges{4}(i);
    im = squeeze(smm(rangestmp{:}));
    if use_spm
      filenameout = [filename sprintf(['%0' num2str(nchars) 'd'],rangestmp{4}) extname];
      v.fname = filenameout;
      spm_write_vol(v,im);
    else
      filenameout = [filename sprintf(['%0' num2str(nchars) 'd'], ...
					 rangestmp{4})];
      if (ndims(im) > 2)
        % Need to have the 3rd coordinate be a scalar so that the image
        % type can be inferred correctly
        pdims = [1 2 ndims(im)+1 3:ndims(im)];
        im = permute(im,pdims);
      end
      dicomwrite(im, filenameout);
    end
    pb_options.progress = i;
  end
  pb_options = progress_bar(pb_options);
  cd(d_current);   % restore old directory
  