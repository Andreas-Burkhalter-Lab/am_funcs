function msprof_nnmf(varargin)
% msprof_nnmf: decompose a set of mass spec runs into m/z peaks
%
% Syntax:
%    msprof_nnmf  % runs with default file names
%    msprof_nnmf(basename)   % uses custom file names
%    msprof_nnmf(options)    % set custom parameter values
%    msprof_nnmf(basename,options)
% 
% This is the main "get started" function in processing mass spec profile
% data. It performs m/z registration, temporal registration (although this
% is complex, and this only does a global shift), and then extracts the m/z
% peaks and calculates their m/z and temporal profiles.  It creates
% files in the current directory that have the following default names:
%    msprof_ecs_filelist: holds the names of the files to process
%    msprof_ecs_registration: holds registration information
%    msprof_ecs_nnmf: results from nonnegative matrix factorization (the
%      mz/ and temporal profiles)
%    msprof_ecs_compounds: a preliminary estimate of the elution-time
%      profiles that correspond to single compounds. "Manual" editing,
%      using msprof_choose_peak_timecourses_gui, is highly recommended.
%
% In general, this function needs to be run twice: once to do registration,
% and once to perform NNMF on the registered data. In both cases, the
% syntax can be as simple as
%   msprof_extract_compounds_set
% On the first run, you will be prompted to define a set of mzXML files to
% process, and optionally (recommended) define a "standard" run to serve
% as a mass-accuracy check and/or calibration as well as a definition of
% the "background." You will also be given the opportunity to pick a set of
% peaks that are reasonably common across your sample set; these will be
% used to correct any m/z drift over time (the registration step). If this
% is not necessary, you can just close the windows to proceed to the last
% step of the first run, which is to define a temporal window of interest.
%
% On the second run, the files will be re-loaded and then processed
% automatically to represent individual peaks as a product of an m/z
% profile and a temporal elution profile, via a sequential single-component
% nonnegative matrix factorization. Intensity "bins" that are
% well-explained by the NNMF are then suppressed from the data for the next
% component. The common isotopologues are also extracted, and an estimate
% of the abundance of the most common elements that have isotopologues
% (e.g., C, O, and S) is made; the peaks corresponding to these
% isotopologues are also suppressed before the next peak is extracted. The
% algorithm starts with the largest peak and works towards smaller ones.
% Finally, after the NNMF, an attempt at identifying peaks corresponding to
% individual compounds in the temporal elution profiles is made. At
% present, this last process should be best viewed as an attempt at a "good
% start" and is unlikely to be satisfactory without manual
% intervention via msprof_choose_peak_timecourses_gui. You can run this GUI
% immediately after completion of this function to see what you've got.
%
% More details about the algorithm can be found in the document msprof.pdf
% in this same directory.
%
% There are numerous options that you can set using the syntax
%    msprof_nnmf(options)
% where options is a structure. Here are the most important fields:
%   isotopeoptions: controls which elements are considered as candidates
%     for the composition of individual m/z peaks. See "isotopes" for a
%     complete description. Your choice here affects how the isotopologue
%     analysis works, because only isotopologues of the reasonably-abundant
%     isotopes among the elements in the "list" are considered. The
%     defaults are described in isotopes. I suspect you must have both
%     'C' and 'S' present, or you will get errors.
%   max_peaks (default Inf): maximum number of m/z peaks that will be
%     extracted.
%   mzregister (default false): if true, the m/z values are first
%     "registered" to compensate for instrument calibration errors. It
%     helps to have a standard compound of known m/z when doing this.
%   nu (default 1/3): power by which to "scale" the intensities so that the
%     errors of the scaled intensities are approximately Gaussian. The
%     errors of the original (unscaled) intensities would be I^(1-nu). The
%     easiest way to calibrate this value is to compare Icalc and Imeasured
%     for the biggest base peak in MSPROF_CHARGE_ISOTOPES_LINEAR. The LTQ
%     seems to fit nu = 1/3 (see also Li et al, Rapid. Comm. Mass Spec.
%      20: 1551–1557, 2006).
%   zscore (default 3): determines threshold via the 'z-score' needed for
%     significant peaks; accepted peaks must be at least z absolute
%     deviations above the median (in counts^nu).
%   peakwidth_mzI (default 5): number of discrete m/z units that define the
%     half-width of a single compound. This should be defined as the
%     "inner" width, i.e., the region that contains most of the power of
%     the signal
%   blankwidth_mzI (default 100) and other options described in
%     msprof_charge_isotopes_linear: FTMS instruments yield peaks with
%     "sidelobes" that can extend well beyond the width of the main peak.
%     This parameter (in conjunction with others as described in
%     msprof_charge_isotopes_linear) defines the extent of the
%     region that, after identifying a peak, will be "blanked" if the
%     NNMF model sufficiently matches the data. Thus, this should be set to
%     the "outer" width of your peaks (i.e., the region containing any
%     observable elevation above baseline). If you want to consider changing
%     this setting, it is highly recommended to look at the ion abundance
%     using logarithmic scaling.
%  files_filename, registration_filename, nnmf_results_filename, and
%    compounds_filename: use these to change the default filenames
%    (described above) for the output.
%
% See also: msprof_choose_peak_timecourses_gui, msprof_charge_isotopes_linear, isotopes.

% Copyright 2009 by Timothy E. Holy

  options = struct;
  basename = 'msprof_ecs';
  for i = 1:length(varargin)
    if isstruct(varargin{i})
      options = varargin{i};
    elseif ischar(varargin{i})
      basename = varargin{i};
    end
  end
  options = default(options,...
    'files_filename',[basename '_filelist.mat'],...
    'mzregistration_filename',[basename '_mzregistration.mat'],...
    'nnmf_results_filename',[basename '_nnmf.mat'],...
    'mzregister',false,...
    'mzregister_thresh',0.1,...
    'mzregister_maxpeaks',10,...
    'peakwidth_mzi',5,...
    'blankwidth_mzi',100,...
    'lookwidth_mzi',10,...
    'nu',1/3,...
    'max_peaks',Inf,...
    'zscore',3,...
    'isotopeoptions',struct,...
    'usemax',true);
  
  % Load the file names
  f = load(options.files_filename);
  n_files = length(f.file);
  file = f.file;
  isotopeinfo = f.isotopeinfo;
  
  %% If necessary, do the m/z registration
  s = struct;
  if options.mzregister
    if (n_files > 1 && ~exist(options.mzregistration_filename,'file'))
      % Load the data at full resolution
      [Ic,~,i2mz] = msprof_load_set(file);
      mzi_max = max(cellfun(@(p) max([p.rowi_max]),Ic));
%       mz = i2mz(1:mzi_max)';
%       col = zeros(n_files,3);  % line color
      % Compute sum across time for each file
      Imz = zeros(mzi_max,n_files);
      for k = 1:n_files
        tmp = ssparse_sum(Ic{k},2);
        Imz(1:length(tmp),k) = tmp;
%         col(k,:) = unique_color(k,n_files);
      end
      % Compute the max intensity across files
      Imz_max = max(Imz,[],2);
      % Find enough peaks that we can unambiguously compare files
      rng0 = [-1 1]*options.lookwidth_mzi;
      rng0 = rng0(1):rng0(2);
      crossterms = zeros(n_files,n_files);  % will be used to determine whether we have enough overlap
      peakList = [];
      isconnected = false;
      while (~isconnected && length(peakList) < options.mzregister_maxpeaks)
        % Pick the biggest peak remaining
        [~,maxIndex] = max(Imz_max);
        rng = rng0+maxIndex;
        skipFlag = false;
        if (any(rng < 1) || any(rng > mzi_max))
          skipFlag = true;
        end
        % Blank it so we don't pick this one again
        rng = rng(rng > 0 & rng <= mzi_max);
        Imz_max(rng) = 0;
        if skipFlag
          continue
        end
        % Store this as a good peak
        peakList(end+1) = maxIndex; %#ok<AGROW>
        % Compute the relative intensity in each file
        thisI = sum(Imz(rng,:),1);
        thisI = thisI/max(thisI);
        % Add to our collection of overlaps
        crossterms = crossterms + thisI'*thisI;
        isconnected = length(connected_components(crossterms >= options.mzregister_thresh)) == 1;
      end
      if ~isconnected
        error('Cannot find overlaps among peaks. Consider skipping registration.');
      end
      % Mutually register these peaks by maximizing overlap
      rng0 = -options.peakwidth_mzi:options.peakwidth_mzi;
      func = @(deltamzi) pairwise_correlation(Imz,peakList,rng0,deltamzi);
      ops = optimset(@fmincon);
      ops.DiffMinChange = 0.1;
      ops.DiffMaxChange = 1;
      meanShiftMzI = fmincon(func,zeros(1,n_files),[],[],ones(1,n_files),0,[],[],[],ops);
      % Plot the results
      fprintf('Registration performed using the following peaks:\n');
      fprintf('%g ',i2mz(peakList));
      fprintf('\n');
      Isnip = snip_interp(Imz,peakList,rng0,meanShiftMzI);
      figure; subplot(1,2,1); plot(Isnip); subplot(1,2,2); plot(meanShiftMzI); xlabel('File #'); ylabel('mzi shift');
      % Save the results
      save(options.mzregistration_filename,'file','meanShiftMzI');
      s.meanShiftMzI = meanShiftMzI;
    else
      s = load(options.mzregistration_filename);
      indx = findainb(file,s.file);
      s.meanShiftMzI = s.meanShiftMzI(indx);
    end
  end
  s = default(s,'meanShiftMzI',zeros(1,n_files));

  %% Perform NNMF to extract m/z and temporal profiles
  %
  % Load registered data
  %
  mzShift = -s.meanShiftMzI;
  mzShift = mzShift - min(mzShift);
  [Ic,mz2i,i2mz,tc] = msprof_load_set(file,struct('mzi_shift',mzShift));
  % Aggregate the data
  I = [Ic{1:n_files}];
  n_scans = length(I);
  n_mzi = max([I.rowi_max]);
  mz = i2mz(1:n_mzi)';
  Imz = zeros(n_mzi,n_files);
  for k = 1:n_files
    if options.usemax
      tmp = ssparse_max(Ic{k},2);
      Imz(1:length(tmp),k) = max(tmp,0);
    else
      tmp = ssparse_sum(Ic{k},2);
      Imz(1:length(tmp),k) = tmp;
    end
  end
  Ip = ssparse_swap(I,n_mzi); % Generate a version quickly indexed by m/z
  %
  % If the user provided a standard formula, determine the shift required
  % bring alignment to the exact mass
  %
  standardMzIShift = 0;
  if ~isempty(f.standardfileindex) && ~isempty(f.standardformula)
    stdImz = ssparse_sum(Ic{f.standardfileindex},2);
    % We require that the standard compound be the largest peak in the
    % sample
    [~,maxIndex] = max(stdImz);
    rng = (-options.lookwidth_mzi:options.lookwidth_mzi) + maxIndex;
    rng = rng(rng >= 1 & rng <= length(stdImz));
    meanMz = sum(stdImz(rng).*mz(rng)) / sum(stdImz(rng));
    standardmass = sum(f.standardformula .* [isotopeinfo.base_mass]);
    mzi = mz2i([meanMz standardmass]); % indices corresponding to masses
    standardMzIShift = diff(mzi); % the index shift required for alignment
  end
  %
  % Set the threshold for peak detection
  %
  Imz_max = max(Imz,[],2);
  if isempty(f.blankindex)
    % Set threshold from main samples
    scaled_max = Imz_max.^options.nu;
  else
    % Set threshold in terms of blank files
    scaled_max = max(Imz(:,s.blankindex),[],2).^options.nu;
  end
  m = median(scaled_max);
  a = mean(abs(scaled_max - m));
  thresh = (m+options.zscore*a)^(1/options.nu);
  save(options.nnmf_results_filename,'file','tc','mz2i','i2mz','standardMzIShift','isotopeinfo','Imz_max','thresh','options')
  %
  % Prepare for iteration over peaks
  %
  % Prepare to split consolidated data by getting breaks in time
  tlen = cellfun(@length,tc);
  breaks_t = [0 cumsum(tlen)];
  % Prepare storage for results
  nnmfresults = struct('mzMeanCorrected',[],'totIntensity',[],'mzCorrected',[],'mzProf',[],...
    'iProf',{cell(1,n_files)},...
    'nnmf_err',[],'Imzrange',[],...
    'charge',[],'chargeI',[],...
    'nI',[],'nIcov',[]);
  rng0 = -options.peakwidth_mzi:options.peakwidth_mzi;
  i = 0;  % this will hold the # of peaks found so far
  fprintf('# processed so far (& fold-threshold of current peak): ');
  % From here on out, Imz becomes an across-files variable and Imz_max
  % becomes a scalar
  if options.usemax
    Imz = Imz_max;
  else
    Imz = sum(Imz(:,1:n_files),2);
  end
  [Imz_max,maxIndex] = max(Imz);
  tic
  tlast = 0;  % keep track of time since last update
  %
  % Iteratively loop over peaks from largest down to threshold, calculating
  % the profile, fitting isotopologues, and then blanking so that we move
  % on to the next peak.
  %
  while (Imz_max > thresh && i < options.max_peaks)
    i = i+1;
    if (i >= length(nnmfresults))
      % Expand the size of the storage by doubling it, so that it doesn't
      % get re-allocated frequently
      nnmfresults(2*i) = nnmfresults(i);
    end
    % Calcualte the range around the base peak, and make sure it fits
    % within the boundaries of the recorded m/z
    rng = maxIndex + rng0;
    rng = rng(rng>=1);
    rng = rng(rng<=length(Imz));
    nnmfresults(i).Imzrange = rng([1 end]);
    % Perform single-component NNMF on the data on the base peak
    Isnip = Ip(rng);
    % Scale so we get approximately Gaussian errors
    for ii = 1:length(Isnip)
      Isnip(ii).value = Isnip(ii).value.^options.nu;
    end
    if options.usemax
      Imz0 = ssparse_sum(Isnip,2);
    else
      Imz0 = Imz(rng).^options.nu;
    end
    [mzProf,iProf,nnmfresults(i).nnmf_err] = ssparse_nnmf1(Isnip,struct('w0',Imz0));
    % Undo the scaling
    mzProf = mzProf.^(1/options.nu);
    iProf = iProf.^(1/options.nu);
    % Massage for output
    iProf(end+1:n_scans) = 0;  % make sure it has full length
    mzProf(end+1:length(rng)) = 0;
    % Report per-"pixel" average error rather than summed error
    nnmfresults(i).nnmf_err = nnmfresults(i).nnmf_err/(n_scans*length(rng));
    % Swap the normalization
    nrm = sqrt(sum(mzProf.^2));
    mzProf = mzProf/nrm;
    nnmfresults(i).mzProf = mzProf;
    iProf = iProf*nrm;
    mzMean = sum(mzProf .* mz(rng))/sum(mzProf);
    indexBase = mz2i(mzMean);
    nnmfresults(i).mzMeanCorrected = i2mz(indexBase + standardMzIShift);
    nnmfresults(i).mzCorrected = i2mz(rng + standardMzIShift);
    nnmfresults(i).totIntensity = sum(iProf)*sum(mzProf);
    for fileIndex = 1:n_files
      rngfile = breaks_t(fileIndex)+1:breaks_t(fileIndex+1);
      nnmfresults(i).iProf{fileIndex} = iProf(rngfile);
    end
    % Do what we can to identify the molecular formula: get the charge
    % and estimate the number of each element using isotopologue
    % abundance. Also blank the base and all isotopologue peaks.
    MpOld = Ip;  %#ok<NASGU> % for debugging: allow you to repeat the previous call
    [nnmfresults(i).charge,nnmfresults(i).chargeI,nnmfresults(i).nI,nnmfresults(i).nIcov,Ip,Imz_changed,index_changed] = msprof_charge_isotopes_linear(mzMean,isotopeinfo,mz2i,Ip,iProf,options);
    if isempty(nnmfresults(i).nI)
      warning('msprof:extract','Aborting because fitting didn''t work')
      break
    end
    % Error checking code:
    if any(Imz(index_changed)*(1+100*eps) < Imz_changed)
      warning('msprof:extract','Something weird happened');
    end
    Ipmx = ssparse_max(Ip(index_changed),2);
    if ~isequal(Imz_changed,Ipmx) && ~(isempty(Imz_changed) && isempty(Ipmx))
      warning('msprof:extract','Mismatch between Imz_changed and the Ip output');
    end
    Imz(index_changed) = Ipmx;
    %
    % Prepare for the next iteration
    %
    Imz_max_last = Imz_max;
    [Imz_max,maxIndex] = max(Imz);
    if (Imz_max > Imz_max_last*(1+100*eps))
      warning('msprof:extract','Something weird happened');
    end
    tnext = toc;
    if (tnext > tlast+10)
      fprintf('%d (%g)...',i,Imz_max/thresh);
      tlast = tnext;
    end
  end
  fprintf('\n')
  toc
  % Truncate any empty space
  nnmfresults = nnmfresults(1:i); %#ok<NASGU>
  randKey = rand;  %#ok<NASGU> % use this to make sure all files go together
  save(options.nnmf_results_filename,'nnmfresults','randKey','-append')
end % msprof_nnmf
  

function [Isnip,rng] = snip_interp(Imz,peakList,rng0,deltamzi)
  % Evaluate Imz around a set of peaks (or list of any m/z), shifting each
  % by an offset deltamzi
  % Imz is n_mz-by-n_files
  % peaksList is a vector
  % rng0 gives the span around each peak (e.g., -5:5)
  % deltamzi is 1-by-n_files, the shift for each file
  n_peaks = length(peakList);
  [n_mz,n_files] = size(Imz);
  % Make the snipped-out range
  rngc = cell(n_peaks,1);
  for peakIndex = 1:n_peaks
    rngc{peakIndex} = rng0(:) + peakList(peakIndex);
  end
  rng = sort(cat(1,rngc{:}));
  % Do the interpolation
  mzi = 1:n_mz;
  Isnip = zeros(length(rng),n_files);
  for fileIndex = 1:n_files
    Isnip(:,fileIndex) = interp1(mzi,Imz(:,fileIndex),rng+deltamzi(fileIndex));
  end
end

function Q = pairwise_correlation(Imz,peakList,rng0,deltamzi)
  Isnip = snip_interp(Imz,peakList,rng0,deltamzi);
  Q = 0;
  n_files = size(Imz,2);
  for i = 1:n_files
    for j = i+1:n_files
      Q = Q + sum(Isnip(:,i).*Isnip(:,j));
    end
  end
  Inorm = sum(Isnip(:).^2);
  Q = -Q/Inorm;
end
