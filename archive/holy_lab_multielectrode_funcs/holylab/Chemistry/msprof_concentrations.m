function s = msprof_concentrations(varargin)
% msprof_concentrations: measure compound concentrations from NNMF and timeSpans
% Syntax:
%   s = msprof_concentrations
%   s = msprof_concentrations(basename)
%   s = msprof_concentrations(options)
%
% The default values for the filenames are:
%   filenameNNMF = 'msprof_ecs_nnmf.mat';
%   filenameCompounds = 'msprof_ecs_compounds.mat';
% On output, s has the following fields:
%   mz: a 1-by-n_compounds vector giving the m/z associated with the
%     compound
%   xrange: a 2-by-n_compounds matrix, time where each column is the
%     beginning/ending coordinate for its elution (if time, in minutes; can
%     also be in %B if standardizing to a gradient)
%   I: is a n_samples-by-n_compounds matrix, giving the total ion current
%     associated with each compound in each sample.
%   xaxisstr: either 'Time (min)' or '%B' depending on whether elution is
%     measured  in time or gradient state, respectively.

% Copyright 2009-2010 by Timothy E. Holy

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
    'nnmf_results_filename',[basename '_nnmf.mat'],...
    'compounds_filename',[basename '_compounds.mat']);
  
  % Choose m/z values that contain timeSpan data
  sc = load(options.compounds_filename);
  sn = load(options.nnmf_results_filename);
  if ~isequal(sc.randKey,sn.randKey)
    error('Keys do not match');
  end
  l = cellfun(@(p) size(p,2),sc.timeSpans);
  keepFlag = l > 0 & ~cellfun(@iscell,sc.timeSpans) & ~cellfun(@isempty,sc.timeSpans);
  tS = sc.timeSpans(keepFlag);
  nnmfresults = sn.nnmfresults(keepFlag);
  mzIndex = find(keepFlag);
  if isfield(sc,'xaxisstr') && isfield(sn,'xaxisstr')
    if ~isequal(sc.xaxisstr,sn.xaxisstr)
      error('x-axes do not match');
    end
    if strcmp(sc.xaxisstr,'Time (min)')
      xc = sn.tc;
    else
      xc = sn.xc;
    end
    s.xaxisstr = sc.xaxisstr;
  else
    s.xaxisstr = 'Time (min)';
    xc = sn.tc;
  end
  
  % Calculate total intensity over timespans
  n_umz = length(nnmfresults);
  n_compounds = sum(l(keepFlag));
  n_files = length(sn.file);
  mz = zeros(1,n_compounds);
  xrange = zeros(2,n_compounds);
  I = zeros(n_files,n_compounds);
  isotopeinfo = isotopes;
  nIsotopes = length(isotopeinfo);
  nI = zeros(nIsotopes,n_compounds);
  curIndex = 1;
  for umzIndex = 1:n_umz
    for cIndex = 1:size(tS{umzIndex},2)
      thisXRange = tS{umzIndex}(:,cIndex);
      mz(curIndex) = nnmfresults(umzIndex).mzMeanCorrected;
      xrange(:,curIndex) = thisXRange;
      iProf = nnmfresults(umzIndex).iProf;
      for fileIndex = 1:n_files
        thisI = iProf{fileIndex};
        thisX = xc{fileIndex};
        flag = thisX >= thisXRange(1) & thisX <= thisXRange(2);
        I(fileIndex,curIndex) = sum(thisI(flag));
      end
      nI(:,curIndex) = nnmfresults(umzIndex).nI(:);
      curIndex = curIndex+1;
    end
  end
  s.mz = mz;
  s.xrange = xrange;
  s.I = I;
  s.nI = nI;
  if exist('msprof_ecs_filelist.mat','file')
    tmp = load('msprof_ecs_filelist.mat','samplenames');
    s.stimtag = tmp.samplenames(1:n_files);
  end
  