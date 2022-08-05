function msprof_extract_compounds_set(varargin)
% msprof_extract_compounds_set: process a set of mass spec profiles
%
% This is the main "get started" function in processing mass spec profile
% data. It performs m/z registration, temporal registration (although this
% is complex, and this only does a global shift), and then extracts the m/z
% peaks and calculates their m/z and temporal profiles.  It creates
% files in the current directory that have the following default names:
%    msprof_ecs_filelist: holds the names of the files to process
%    msprof_ecs_mzregistration: holds registration information
%    msprof_ecs_nnmf: results from nonnegative matrix factorization (the
%      m/z and temporal profiles of peaks)
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
%    msprof_extract_compounds_set(options)
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
%   z (default 5): determines threshold via the 'z-score' needed for
%     significant peaks; accepted peaks must be at least z absolute
%     deviations above the median (in log(counts)).
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
%  peak_modulation (default 2): the factor by which peaks must exceed
%    their adjacent troughs in order to be considered to be a well-isolated
%    temporal peak. Setting this lower will give you more
%    automatically-identified single-compound peaks (at the possible risk
%    of being distracted by noise).
%  files_filename, registration_filename, nnmf_results_filename, and
%    compounds_filename: use these to change the default filenames
%    (described above) for the output.
%
% See also: msprof_choose_peak_timecourses_gui, msprof_charge_isotopes_linear, isotopes.

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
    'peak_modulation',2,...
    'isotopeoptions',struct,...
    'files_filename',[basename '_filelist.mat'],...
    'mzregistration_filename',[basename '_mzregistration.mat'],...
    'nnmf_results_filename',[basename '_nnmf.mat'],...
    'compounds_filename',[basename '_compounds.mat']);
  % Options for identifying compounds from the intensity timecourse
  options = default(options,'t_bin',0.1,'t_min',0.5,'adevfac',2,'standardfac',10);
  % Prepare for analysis of isotopologues
  isotopeinfo = isotopes(options.isotopeoptions);
  
  %% Get the files
  if ~exist(options.files_filename,'file')
    file = UIGetFiles('*.mzXML','Select your files');
    if isempty(file)
      return
    end
    samplenames = cell(1,length(file));
    for sampleIndex = 1:length(file)
      [~,samplenames{sampleIndex}] = fileparts(file{sampleIndex});  % default is the file name
    end
    [blankindex,ok] = listdlg('ListString',samplenames,'ListSize',[400 400],'SelectionMode','multiple','PromptString','Choose any of these that are blanks. If none, click cancel'); %#ok<NASGU>
    if ~ok
      blankindex = [];
    end
    standardformula = []; %#ok<NASGU>
    standardcharge = []; %#ok<NASGU>
    standardfileindex = []; %#ok<NASGU>
    buttonname = questdlg('Did you dope any of your samples with a standard?');
    switch buttonname
      case 'Cancel'
        return
      case 'No'
        % do nothing
      case 'Yes'
        buttonname = 'No';
        while isequal(buttonname,'No')
          fstr = inputdlg({'Enter the molecular formula of the standard base ion (e.g., "C21 H29 O7 S")','Enter charge of this ion'},'Enter chemical formula');
          if ~isempty(fstr)
            standardformula = string2formula(fstr{1},isotopeinfo);
            standardcharge = abs(str2double(fstr{2}));
            standardmass = sum(standardformula .* [isotopeinfo.base_mass]);
            buttonname = questdlg(sprintf('The compound you entered has a base peak of %.4f. Is this correct?',standardmass/standardcharge));
          else
            standardformula = [];
            break
          end
        end
        if isequal(buttonname,'Cancel')
          standardformula = [];
        end
        standardfileindex = listdlg('ListString',file,'SelectionMode','single','PromptString','Select a file in which the standard is the dominant peak (if none, hit cancel)');
    end  % getting filenames and standards
    
    % Optional file renaming
    uiwait(msgbox('If you want to "rename" the samples, pick the .csv (text) file with the names (column 1: file names; column 2: new names). Or, hit cancel to use the existing filenames.','Renaming','none','modal'));
    [namefile,pth] = uigetfile('*.csv');
    if (~isempty(namefile) && ~isnumeric(namefile))
      [fid,msg] = fopen([pth filesep namefile],'r');
      if (fid < 0)
        error(msg);
      end
      tline = cell(0,1);
      while (1)
        tline{end+1} = fgets(fid); %#ok<AGROW>
        if (tline{end} == -1)
          break
        end
      end
      fclose(fid);
      n_lines = length(tline)-1;
      namelookup = cell(n_lines,2);
      for lineIndex = 1:n_lines
        mtch =  regexp(tline{lineIndex},'".*?"','match');
        for i = 1:2
          namelookup{lineIndex,i} = mtch{i};
          % Strip initial and final ", if present
          if (namelookup{lineIndex,i}(1) == '"')
            namelookup{lineIndex,i} = namelookup{lineIndex,i}(2:end-1);
          end
        end
      end
      % Replace matching sample names
      for sampleIndex = 1:length(file)
        matchIndex = strmatch(samplenames{sampleIndex},namelookup(:,1),'exact');
        if ~isempty(matchIndex)
          samplenames{sampleIndex} = namelookup{matchIndex,2};
        end
      end
    end
    % Save everything
    save(options.files_filename,'file','samplenames','blankindex','standardfileindex','standardformula','standardcharge','isotopeinfo');
  else
    % Information about files and standards was already present
    load(options.files_filename);
  end

  %% Do the NNMF
  msprof_nnmf(options);
  
  %% Do temporal registration
  snnmf = load(options.nnmf_results_filename);
  gc = set_gradients(file);
  tShift = msprof_temporal_register(snnmf.nnmfresults(1:20),snnmf.tc,gc);
  save([basename '_tregistration'],'tShift','file');
  
  return
  
  %% Find temporal peaks and quantify intensity in them
  % Prepare to resample evenly on the time axis
  n_mz = length(nnmfresults);
  n_rs = ceil(diff(s.trange)/options.t_bin);
  options.n_smooth = 1;
  options.n_min = round(options.t_min / options.t_bin);
  [trs,Mrs] = resample_evenly_preserve_integral(tc,n_rs);
  Irs = zeros(n_rs,length(Mrs));
  % Compute intensity & variability statistics on intervals with peaks in them
  pvar = cell(1,n_mz);
  Itot = cell(1,n_mz);
  ItotBlank = cell(1,n_mz);
  timeSpans = cell(1,n_mz);
  tic;
  tlast = 0;
  fprintf('Processing %d m/z values: ',n_mz);
  for mzIndex = 1:n_mz
    if (toc - tlast > 5)
      fprintf('%d...',mzIndex);
      tlast = toc;
    end
    % Resample evenly
    for fileIndex = 1:n_files
      Irs(:,fileIndex) = Mrs{fileIndex}*nnmfresults(mzIndex).iProf{fileIndex}(:);
    end
    if ~isempty(standardfile)
      Irs(:,end) = Mrs{end}*nnmfresults(mzIndex).iProfBlank(:);
    end
    % Find the temporal peaks, using just the "real" samples
    [rngc,pvar{mzIndex}] = timespan_frompeaks(Irs(:,1:n_files),options);
    % Compute intensity in peaks
    n_spans = length(rngc);
    Itottmp = zeros(n_spans,size(Irs,2));
    for tsIndex = 1:n_spans
      Itottmp(tsIndex,:) = sum(Irs(rngc{tsIndex},:),1);
    end
    Itot{mzIndex} = Itottmp(:,1:n_files);
    if ~isempty(standardfile)
      ItotBlank{mzIndex} = Itottmp(:,end);
    end
    timeSpans{mzIndex} = [trs(cellfun(@(p) p(1),rngc)); trs(cellfun(@(p) p(end),rngc))];
  end
  fprintf('done.\n');
  
  % Determine a threshold on the intensity
  It = cat(1,Itot{:});
  ItBlank = cat(1,ItotBlank{:});
  mIt = max(It,[],2);
  Ifig = figure; hist(log10(mIt+1),200); title('Distribution of peak intensities')
  if ~isempty(ItBlank)
    figure; hist(log10(ItBlank+1),200); title('Distribution of standard (blank) intensities')
    lt = log10(ItBlank(ItBlank > 0));
    mlt = median(lt);
    alt = mean(abs(lt - mlt));
    thresh = 10^(mlt+options.adevfac*alt);
  else
    yl = ylim;
    xl = xlim;
    hline = line([1 1]*mean(xl),yl,'Color','r','LineWidth',2);
    drag_line(hline);
    setappdata(Ifig,'threshH',hline);
    set(Ifig,'CloseRequestFcn',@(obj,event) uiresume_wrapper(obj,event,Ifig));
    msgbox('Drag the line to select the threshold, and close window when done','Select threshold','modal');
    if ishandle(Ifig)
      uiwait(Ifig)
    end
    thresh = get(hline,'XData');
    thresh = 10^thresh(1);
    delete(Ifig);
  end    
  keepFlag = mIt > thresh;
  fprintf('Out of %d identified "compounds", %d exceeded the threshold.\n',length(keepFlag),sum(keepFlag));
  if ~isempty(ItBlank)
    keepFlag = keepFlag & mIt > options.standardfac*ItBlank;
    fprintf('Of these, %d compounds were %g-fold larger than in the standard.\n',sum(keepFlag),options.standardfac);
  end

  % Keep only those peaks that exceed the threshold
  l = cellfun(@(p) size(p,2),timeSpans);
  breaks = [0 cumsum(l)];
  for mzIndex = 1:n_mz
    kF = keepFlag(breaks(mzIndex)+1:breaks(mzIndex+1));
    timeSpans{mzIndex} = timeSpans{mzIndex}(:,kF);
    Itot{mzIndex} = Itot{mzIndex}(kF,:)';
    pvar{mzIndex} = pvar{mzIndex}(kF);
  end
  
  save(options.compounds_filename,'timeSpans','randKey','Itot','pvar','thresh');
  
end

function mspecs_choosepeaks(hax,~,options)
  hfig = get_parent_fig(hax);
  cp = get(hax,'CurrentPoint');
  col = getappdata(hfig,'col');
  Mcp = getappdata(hfig,'Mcp');
  tc = getappdata(hfig,'tc');
  Mmz = getappdata(hfig,'Mmz');
  mz = getappdata(hfig,'mz');
  [n_mz,n_files] = size(Mmz);
  [~,~,vertexIndex] = findpoint(cp(1,1:2),hax);
  rng = vertexIndex + [-1,1]*options.lookwidth_mzI;
  rng(1) = max(1,rng(1));
  rng(2) = min(n_mz,rng(2));
  rng = rng(1):rng(2);
  hfignew = figure;
  subplot(1,2,1)
  for k = 1:n_files
    line(mz(rng),Mmz(rng,k),'Color',col(k,:));
  end
  xlabel('m/z')
  subplot(1,2,2)
  tpk = cell(1,length(Mcp));
  for k = 1:length(Mcp)
    tpk{k} = ssparse_sum(Mcp{k}(rng),1);
    tpk{k}(end+1:length(tc{k})) = 0;
    line(tc{k},tpk{k},'Color',col(k,:));
  end
  xlabel('Time (min)')
  setappdata(hfignew,'mainFigure',hfig);
  setappdata(hfignew,'timeCourse',tpk);
  setappdata(hfignew,'rng',rng);
  peakfigH = getappdata(hfig,'peakfigH');
  peakfigH(end+1) = hfignew;
  setappdata(hfig,'peakfigH',peakfigH);
end
  
function mspecs_register(hfig,~,options)
  peakfigH = unique(getappdata(hfig,'peakfigH'));
  peakfigH = peakfigH(ishandle(peakfigH));
  if ~isempty(peakfigH)
    n_figs = length(peakfigH);
    % Do the m/z shift
    Mmz = getappdata(hfig,'Mmz');
    n_files = size(Mmz,2);
    xcm = zeros(n_figs,n_files);
    w = zeros(n_figs,n_files);
    for k = 1:n_figs
      rng = getappdata(peakfigH(k),'rng');
      thisMmz = Mmz(rng,:);
      figure(peakfigH(k));
      x = repmat((1:size(thisMmz,1))',1,n_files);
      tmp = sum(thisMmz.*x,1) ./ sum(thisMmz,1);
      xcm(k,:) = tmp - mean(tmp);
      w(k,:) = sum(thisMmz,1);  % in computing the averages, weight by intensity
    end
    figure
    plot(xcm')
    xlabel('Filenumber')
    ylabel('m/z index shift')
    meanShiftMzI = sum(xcm.*w,1) ./ sum(w,1); %#ok<NASGU>
    % Do the temporal shift---look for the moment when it first crosses
    % half-max
    tc = getappdata(hfig,'tc');
    tcm = zeros(n_figs,n_files);
    for k = 1:n_figs
      Mt = getappdata(peakfigH(k),'timeCourse');
      for kk = 1:n_files
        mxVal = max(Mt{kk});
        indx = find(Mt{kk} > mxVal/2,1,'first');
        tcm(k,kk) = tc{kk}(indx);
      end
    end
    tcm = tcm - repmat(median(tcm,2),1,n_files);
    figure
    plot(60*tcm')
    xlabel('Filenumber')
    ylabel('Scan time shift (s)')
    meanShiftT = sum(tcm.*w,1)./sum(w,1); %#ok<NASGU>
    save(options.registration_filename,'meanShiftMzI','meanShiftT');
    close(peakfigH);
    delete(hfig)
  else
    delete(hfig)
  end
end

function mspecs_choose_trange(hfig,eventdata,options) %#ok<INUSL>
  hline = getappdata(hfig,'trangeH');
  xd = get(hline,'XData');
  trange = sort([xd{1}(1) xd{2}(1)]); %#ok<NASGU>
  if exist(options.registration_filename,'file')
    save(options.registration_filename,'trange','-append');
  else
    save(options.registration_filename,'trange');
  end
  delete(hfig)
end

function uiresume_wrapper(~,~,figh)
  uiresume(figh);
end

function [rngc,pvar] = timespan_frompeaks(Irs,options)
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
  % Make sure there's a real "peak" over each range
  mx = cellfun(@(r) max(Irssum(r)),rngc);
  left = cellfun(@(r) Irssum(r(1)),rngc);
  right = cellfun(@(r) Irssum(r(end)),rngc);
  keepFlag = mx > options.peak_modulation*max([left;right],[],1);
  rngc = rngc(keepFlag);
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
end