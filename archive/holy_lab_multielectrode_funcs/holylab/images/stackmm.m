classdef stackmm
%  STACKMM:  memory-mapped image stack access class
%  Syntax:
%    smm = stackmm(basename)
%    smm = stackmm(basename, 'bias', b)
%    smm = stackmm(files)
%  where
%    basename is a string containing the name of the image header, and
%    image data as a 4D tensor, smm(xindx, yindx, zindx, tindx), where
%    each of these indices can be ':' as well as particular values.  This
%    allows very fast and flexible slicing of the image data.  Be careful
%    not to ask for more than memory permits!
%  and where
%    files is a cell array containing names (or paths if in different
%    folders than the current working directory) of the image headers, as
%    described above.  If more than one file are to be opened 
%    simultaneously, the filenames should be supplied as a cell array. Note
%    that if acquisition of any individual file was prematurely truncated,
%    the .imagine file should be changed to reflect the actual data range.
%    On a similar note, it is not currently necessary to change the
%    requested frames/stacks value.  Together nframes and frames_per_stack 
%    provide the original value of the requested stacks.  Compatible 
%    headers will be recorded by the same camera as well as having a number 
%    of other header fields in common (all except {'date', 'wholeheader', 
%    'header_filename','idle','nframes','nstacks', 'stacktime',
%    'stim_lookup','stim_labels','trial_lookup'} must be identical).  
% 
%  The 'bias' syntax allows you to supply a value for the camera's bias;
%  this can be either a scalar value (for uniform bias) or a full image
%  (if the bias varies across pixels).  This value will be subtracted
%  from each pixel before returning; you can turn bias subtraction off
%  by
%      smm.bias = [];
%  and back on by
%      smm.bias = b.
%  You can also get additional info from the following:
%    smm.size:  the size of the 4D stack
%    smm.header:  the image file header
%    smm.type:  the raw data type, e.g. 'uint16'
%    smm.filename:  the name of the file on disk
%    smm.bias:  the current bias
%  Users can also assign the following fields (although these data fields 
%  are currently ONLY maintained while the stackmm object is loaded into 
%  the workspace, ie the user must maintain the data offline and reload on
%  subsequent desired use):
%    smm.badframes: a logical(size(3),size(4)) matrix denoting, by logical 
%     true, whether each frame(row vector) in each stack(column vector)
%     is known to contain bad data.  For example, this field may be loaded
%     from the logical isbad matrix output of find_bad_frames.m for
%     subsequent exlusion of "bad" data from registration.  In this
%     example, after registration, image shift (or registration) results in
%     the badframes and badpixel data fields being in accurate with regard
%     to which data they refer.  For this reason badframes is not
%     permanently associated with the header file.  
%    smm.badpixels: a logical (size(1),size(2)) matrix denoting, by logical
%     true, whether a pixel in each frame is known to be bad.  For example,
%     in uint16 data recorded from a camera, a pixel value which is always
%     zero can now be annoted using this field.
%
%  
%
%  Copyright 2006, 2011 by Timothy E. Holy and Gary F. Hammen
    
  properties (GetAccess = public, SetAccess = private)
    header = struct;
    filename = '';
    bias = [];
    badframes = [];
    badpixels = [];
  end
  properties (Dependent = true, SetAccess = private)
    size
    type
  end
  properties (Access = private)
    mm = [];     % the memory-mapped object(s)
    sz3 = [];
    n_stacks = [];
  end
  methods
    %% Constructor
    function smm = stackmm(basename, varargin)
      if isa(basename, 'stackmm')
        % Copy Constructor
        smm = basename;
      else
        % Create a new stackmm
        smm.filename = basename;
        if ischar(basename)
          basename = {basename};
        end
        n_files = length(basename);
        smm.n_stacks = zeros(1,n_files);
        for fileIndex = n_files:-1:1
          h(fileIndex) = imreadheader(basename{fileIndex});
        end
        % Check to see that the headers are consistent
        fields_changeable = {'date', 'wholeheader', 'header_filename','idle','nframes','nstacks','stacktime','stim_lookup','stim_labels','trial_lookup'};
        fields_changeable = intersect(fields_changeable,fieldnames(h));
        h_fixed = rmfield(h,fields_changeable);
        issame = arrayfun(@(x) isequal(x,h_fixed(1)), h_fixed);
        if any(~issame)
          error('One or more non-changeable fields differs from one file to the next');
        end
        % Create a merged header
        h1 = h(1);
        if (n_files > 1)
          % For the changeable fields, many will be concatenations of the
          % values for each file
          h1.date = {h.date};
          h1.wholeheader = {h.wholeheader};
          h1.header_filename = {h.header_filename};
          h1.idle = [h.idle];
          h1.nframes = [h.nframes];
          h1.nstacks = [h.nstacks];
          % Convert stacktime to an absolute time relative to the start of
          % the experiment
          tstart = zeros(1,n_files);
          for fileIndex = 1:n_files
            [s1,s2] = strtok(h(fileIndex).date,'T');
            tstart(fileIndex) = datenum([s1 ' ' s2(2:end)]);
            h(fileIndex).stacktime = h(fileIndex).stacktime + (24*3600)*(tstart(fileIndex) - tstart(1));
          end
          h1.stacktime = [h.stacktime];
          % Handle the stimuli. If they are all the same, then preserve them
          % as they are
          if ~all(cellfun(@isempty,{h.stim_lookup}))
            eq = true;
            h1.stim_lookup = [];
            for fileIndex = 2:n_files
              eq = isequal(h(fileIndex).stim_labels,h(1).stim_labels);
              if ~eq
                break
              end
            end
            if eq
              h1.stim_labels = h(1).stim_labels;
              for fileIndex = 1:n_files
                tmp = h(fileIndex).stim_lookup(1:h(fileIndex).nstacks);
                h1.stim_lookup = [h1.stim_lookup(:); tmp(:)];
              end
            else
              [h1.stim_labels] = unique([h.stim_labels]);
              for fileIndex = 1:n_files
                map = [0 findainb(h(fileIndex).stim_labels,h1.stim_labels)];
                h(fileIndex).stim_lookup = map(h(fileIndex).stim_lookup+1);
              end
              [h1.stim_lookup] = [h.stim_lookup];
            end
            if isfield(h1,'trial_lookup')
              h1 = rmfield(h1,'trial_lookup');
            end
          end
        end
        if ~isfield(h1,'pixel_spacing')
          if isfield(h1,'um_per_pixel_xy')
            if isfield(h1,'piezo_start_stop')
              h1.pixel_spacing = [[1 1]*h1.um_per_pixel_xy ...
                diff(h1.piezo_start_stop)/(h1.frames_per_stack-1)];
            else
              h1.pixel_spacing = [1 1]*h1.um_per_pixel_xy;
            end
          end
        end
        % Store the header
        smm.header = h1;
        smm.sz3 = [h1.width h1.height h1.frames_per_stack];
        % Handle the bias
        if ~isempty(varargin)
          for i = 1:2:length(varargin)
            switch varargin{i}
              case 'bias'
                smm.bias = varargin{i+1};
              otherwise
                error(['Parameter ' varargin{i} ' not recognized'])
            end
          end
        end
        % Set up the memory-mapped file(s). Note that because one can't memory
        % map more than 2G worth of memory at any time, the approach taken
        % here is to map a single stack at a time; the offset is then moved
        % to map particular stacks
        smm.mm = cell(1,n_files);
        n_bytes_per_stack = prod(smm.sz3)*sizeof(smm.header.prec);
        for fileIndex = 1:n_files
          camfilename = stackheaderfilename2camfilename(basename{fileIndex},smm.header);
          if exist(camfilename,'file')
            % Check the file size
            dirtmp = dir(camfilename);
            nstacks = floor(dirtmp.bytes/n_bytes_per_stack);
            if (nstacks < smm.header.nstacks(fileIndex))
              warning('stackmm:smallfile','File %s smaller than expected, changing the number of stacks to %d',camfilename,nstacks);
              smm.header.nstacks(fileIndex) = nstacks;
              h(fileIndex).stacktime = h(fileIndex).stacktime(1:nstacks);
            end
            % Memory-map the file
            smm.mm{fileIndex} = memmapfile(camfilename,...
              'format',{smm.type,smm.sz3,'im'},...
              'repeat',1);
          else
            smm.mm{fileIndex} = [];
            warning('stackmm:filenotfound','.cam file not found');
          end
        end
        smm.header.stacktime = [h.stacktime];
      end
    end
    %% Get methods using the class interface
    function typeout = get.type(smm)
      typeout = smm.header.prec;
    end
    function sz = get.size(smm)
      sz = [smm.sz3 sum(smm.header.nstacks)];
      if strcmp(smm.header.camera,'DV8285_BV')
        % The OCPI1 Andor camera has a line of bad pixels; we'll discard
        % it, and so here we "lie" about the image size to compensate
        sz(1) = sz(1)-1;
      end
    end
    function im = subsref(smm,s)
      switch s.type
        case '()'
          %% Get image data
          % We have to parse indices like ':',
          % check bounds, and then loop over the stack number (4th index)
          % because of the limitation that we memory-map only a single stack at
          % a time.
          sz = [smm.sz3 sum(smm.header.nstacks)];
          stacksize = sizeof(smm.type)*prod(smm.sz3);
          c = cell(1,4);
          l = zeros(1,4);
          % When a single .cam files has been generated from multiple
          % files, the compound flag will side-step the change in mm files
          % that would otherwise be indicated by crossing file boundaries in
          % cumstacks number.
          compound = 0;
          if isfield(smm.header, 'compound_cam')
            compound = smm.header.compound_cam;
          end
          for i = 1:4
            if strcmp(s.subs{i},':')
              c{i} = 1:smm.size(i);
              l(i) = smm.size(i);
            else
              c{i} = s.subs{i};
              l(i) = length(c{i});
              if any(c{i} > sz(i))
                error('ImageStack:outofrange',...
                  'Stack dimensions are only [%d %d %d %d]',sz(1),sz(2),sz(3),sz(4));
              end
            end
          end
          % The OCPI1 Andor camera has a line of bad pixels, discard
          if strcmp(smm.header.camera, 'DV8285_BV')
            if length(c{1}) > 1003
              c{1} = 1:1003;
              l(1) = length(c{1});
            end
          end
          im = zeros(l,smm.type);
          cumstacks = cumsum(smm.header.nstacks); % the cumulative # of stacks/file
          % Loop over the stack index
          for i = 1:length(c{4})
            % Determine which file we should be looking at
            cstack_tot = c{4}(i);
            if true(compound)
              cfile = 1;
              cstack = cstack_tot;
            else
              cfile = find(cstack_tot <= cumstacks,1,'first');
              if (cfile > 1)
                cstack = cstack_tot - cumstacks(cfile-1);
              else
                cstack = cstack_tot;
              end
            end
            % Adjust the offset to move to the current stack
            smm.mm{cfile}.offset = stacksize*(cstack-1);
            % Load the data
            imtmp = smm.mm{cfile}.data.im(c{1:3});
            if ~isempty(smm.bias)
              if isscalar(smm.bias)
                imtmp = imtmp - smm.bias;
              else
                for j = 1:length(c{3});
                  imtmp(:,:,j) = imtmp(:,:,j) - smm.bias(c{1:2});
                end
              end
            end
            % Integrate this stack into the total output
            im(:,:,:,i) = imtmp;
          end
        case '.'
          %% Other "get" methods that use subsref
          % These are necessary due to a limitation in Matlab's class
          % implementation; if you supply a subsref method so you can
          % overload parentheses (which we do for getting image data), then
          % you are also required to handle all other possible subsref
          % calls on your own.
          switch s.subs
            case 'header'
              im = smm.header;
            case 'size'
              im = smm.size;
            case 'type'
              im = smm.type;
            case 'filename'
              im = smm.filename;
            case 'bias'
              im = smm.bias;
            case 'badframes'
              im = smm.badframes;
            case 'badpixels'
              im = smm.badpixels;
            otherwise
              error(['Field ' s.subs ' not recognized']);
          end
      end
    end
    %% Assignment methods
    % These support the user setting "badframes" and "badpixels" after
    % construction
    function obj = subsasgn(obj,s,b)
      switch s.type
        case '.'
          switch s.subs
            case 'badframes'
              if ~isequal(size(b),[obj.sz3(3) sum(obj.header.nstacks)])
                error('The dimensions of the badframes matrix does not match');
              end
              obj.badframes = logical(b);
            case 'badpixels'
              if ~isequal(size(b),obj.size(1:2))
                error('The dimensions of the badpixels matrix does not match');
              end
              obj.badpixels = logical(b);
            otherwise
              error(['Field ' s.subs ' not recognized']);
          end
        otherwise
          error(['subsasgn: ' s.type ' not recognized']);
      end
    end
  end
end