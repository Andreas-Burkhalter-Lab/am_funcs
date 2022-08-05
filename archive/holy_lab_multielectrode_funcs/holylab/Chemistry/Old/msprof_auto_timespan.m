function msprof_auto_timespan(options)
  if (nargin < 1)
    options = struct;
  end
  options = default(options,'t_bin',0.1,'t_min',0.5,'adevfac',3,'filein','msprof_ecs_nnmf.mat','fileout','msprof_ecs_compounds.mat','filelist','msprof_ecs_filelist.mat','blankfac',10);
  
  inputdata = load(options.filein);
  n_mz = length(inputdata.totIntensity);
  n_files = length(inputdata.file);
  tlen = cellfun(@length,inputdata.tc);
  tbreaks = [0 cumsum(tlen)];
  tmin = cellfun(@(p) p(1),inputdata.tc);
  tmax = cellfun(@(p) p(end),inputdata.tc);
  n_rs = ceil(max(tmax)-min(tmin))/options.t_bin;
  options.n_smooth = 1;
  options.n_min = round(options.t_min / options.t_bin);

  [trs,Mrs] = resample_preserve_integral(inputdata.tc,n_rs);
  Irs = zeros(n_rs,n_files);

  % Compute intensity & variability statistics on intervals with peaks in them
  pvar = cell(1,n_mz);
  Itot = cell(1,n_mz);
  timeSpans = cell(1,n_mz);
  tic;
  tlast = 0;
  fprintf('Processing %d m/z values: ',n_mz);
  for mzIndex = 1:n_mz
    if (toc - tlast > 5)
      fprintf('%d...',mzIndex);
      tlast = toc;
    end
    for fileIndex = 1:n_files
      thisI = inputdata.iProf{mzIndex}(tbreaks(fileIndex)+1:tbreaks(fileIndex+1));
      Irs(:,fileIndex) = Mrs{fileIndex}*thisI;
    end
    [rngc,pvar{mzIndex}] = mat_frompeaks(Irs,options);
    n_spans = length(rngc);
    Itottmp = zeros(n_spans,n_files);
    for tsIndex = 1:n_spans
      Itottmp(tsIndex,:) = sum(Irs(rngc{tsIndex},:),1);
    end
    Itot{mzIndex} = Itottmp;
    timeSpans{mzIndex} = [trs(cellfun(@(p) p(1),rngc)); trs(cellfun(@(p) p(end),rngc))];
  end
  fprintf('done.\n');
  
  % Assign a threshold on the intensity
  It = cat(1,Itot{:});
  mIt = max(It,[],2);
  figure; hist(log10(mIt+1),200); title('Distribution of peak intensities')
  lmIt = log10(mIt(mIt > 0));
  mlmIt = median(lmIt);
  almIt = mean(abs(lmIt - mlmIt));
  thresh = 10^(mlmIt+options.adevfac*almIt);
  keepFlag = mIt > thresh;
  fprintf('Out of %d identified "compounds", %d exceeded the threshold',length(keepFlag),sum(keepFlag));
  
  % Keep only those peaks that exceed the threshold
  l = cellfun(@(p) size(p,2),timeSpans);
  breaks = [0 cumsum(l)];
  for mzIndex = 1:n_mz
    kF = keepFlag(breaks(mzIndex)+1:breaks(mzIndex+1));
    timeSpans{mzIndex} = timeSpans{mzIndex}(:,kF);
    Itot{mzIndex} = Itot{mzIndex}(kF,:)';
    pvar{mzIndex} = pvar{mzIndex}(kF);
  end

  % Look for the same peaks in the blank sample
  fl = load(options.filelist);
  if isfield(fl,'blankfile')
    n_killed = 0;
    n_original = 0;
    rf = load('msprof_ecs_registration');
    mzShift = -rf.meanShiftMzI;
    mzShift = mzShift - min(mzShift);
    [Mc,mz,tc] = msload_set([fl.file(1) fl.blankfile],mzShift([1 end]),-rf.meanShiftT([1 end]),rf.trange); % load first to insure same m/z
    if ~isequal(mz,inputdata.mz)
      error('Mismatch in m/z');
    end
    Mp = Mc{2}';
    for mzIndex = 1:n_mz
      if isempty(timeSpans{mzIndex})
        continue
      end
      rngmz = inputdata.breaks(mzIndex)+1:inputdata.breaks(mzIndex+1);
      thisM = Mp(:,rngmz);
      thisI = sum(thisM .* repmat(inputdata.mzProf{mzIndex}(:),1,size(thisM,2)),1);
      kF = false(1,length(timeSpans{mzIndex}));
      for tIndex = 1:length(timeSpans{mzIndex})
        rngt = tc{2} >= timeSpans{mzIndex}(1,tIndex) & ...
          tc{2} <= timeSpans{mzIndex}(2,tIndex);
        Iblank = sum(thisI(rngt));
        kF(tIndex) = any(Itot{mzIndex}(tIndex,:) > Iblank*options.blankfac);
      end
      timeSpans{mzIndex} = timeSpans{mzIndex}(:,kF);
      Itot{mzIndex} = Itot{mzIndex}(kF,:);
      pvar{mzIndex} = pvar{mzIndex}(kF);
      n_killed = n_killed + sum(~kF);
      n_original = n_original + length(kF);
    end
    fprintf('%d time spans killed via comparison to the blank sample',n_killed)
  end

  % Save the results
  randKey = inputdata.randKey;
  save(options.fileout,'timeSpans','randKey','Itot','pvar','thresh');
  

function [rngc,pvar] = mat_frompeaks(Irs,options)
  %% Find peaks and their regions of support in the across-sample summed intensity
  Irss = imfilter_gaussian(Irs,[options.n_smooth 0]);
  Irssum = sum(Irss,2);
  map = imflow(Irssum);
  mapOld = map;
  map = map(map);
  while (~isequal(map,mapOld))
    mapOld = map;
    map = map(map);
  end
  [umap,tmp,index] = unique(map);
  rngc = agglabel(index);
  rngc = rngc(2:end); % truncate edges
  l = cellfun(@length,rngc);
  rngc = rngc(l >= options.n_min); % keep intervals that are "long enough"
  n_intervals = length(rngc);
  if (nargout < 2)
    return
  end

  %% Check each region for "consistency"
  % (the percent of variance explained by the first SVD component, using "registered" data)
  n_files = size(Irs,2);
  pvar = zeros(1,n_intervals);
  for i = 1:n_intervals
    rng = rngc{i};
    % Compute the "center of mass" of each trace within each timespan
    x = repmat((1:length(rng))',1,n_files);
    sI = sum(Irs(rng,:),1);
    sxI = sum(Irs(rng,:).*x,1);
    mx = sxI ./ (sI+1);
    midpoint = length(rng)/2;
    mx(mx == 0) = midpoint;  % for traces with 0 intensity, set center of mass equal to the center of interval
    % "Register" each to its center of mass
    Iresnip = zeros(length(rng),n_files);
    for k = 1:n_files
      resniprng = rng + round(mx(k) - midpoint);
      if resniprng(1) < 1 || resniprng(end) > size(Irs,1)
        resniprng = rng;
      end
      Iresnip(:,k) = Irs(resniprng,k);
    end
    % Do an SVD on Iresnip and see how much is in the first component
    [U,S] = svd(Iresnip,'econ');
    s = diag(S);
    pvar(i) = s(1)^2 / sum(s.^2);
  end
  
