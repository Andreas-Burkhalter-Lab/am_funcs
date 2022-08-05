function [charge,chargeI,n,ncov,Iout,Isum_changed,index_changed] = msprof_charge_isotopes_linear(mzBase,isotopeinfo,mz2i,I,varargin)
% MSPROF_CHARGE_ISOTOPES_LINEAR: determine charge and element composition from isotopologues
%
% This function is valid only for single-substitution isotopologues.  For
% small molecules, these are the most abundant, and thus provide a basis for
% fast fitting because one can use linear algebra.  For more complete
% distributions, see ISOTOPOLOGUES.
%
% This function also supports "blanking," meaning removing/suppressing the
% well-described peaks from the total intensity profile. This is useful in
% contexts in which one wants to successively analyze starting from
% the largest peak and working towards the smaller peaks.
%
% Syntax:
%   [charge,n] = msprof_charge_carbons_sulfurs(mzBase,isotopeinfo,mz2i,I)
%   [charge,n] = msprof_charge_carbons_sulfurs(mzBase,isotopeinfo,mz2i,I,sProf)
%   [charge,n] = msprof_charge_carbons_sulfurs(...,options)
%   [charge,n,ncov,Iout] = msprof_charge_carbons_sulfurs(...)
% where
%   mzBase is the mean m/z value of the base compound;
%   isotopeinfo is data from ISOTOPES;
%   mz2i is a function that converts m/z values into indices
%   I is the intensity "sparse structure matrix," in row-major ordering
%   sProf (optional) is the profile with respect to scan #;
%   options may have the following fields:
%     peakwidth_mzI (default 5): the number of mzi units on either side of the
%       mean m/z to include in the peak
%     blankwidth_mzI (default 100): the number of mzi units over which
%       blanking is to occur (if Iout is requested)
%     blankfactor (default 10) and blankmax (default 1e-3): for any "pixel,"
%       if the measured intensity is lower than
%            blankfactor * (Icalc + blankmax*max(Icalc)),
%       that "pixel" is set to zero in Iout
%     zmax (default 3): the maximum charge to consider
%     dmzi_max (default 3): the maximum number of mzi units allow to slide
%       for "registration"
%     plot (default false): if true, plots temporal data for isotopologues
%     usemax (default false): if true, Isum_changed is the max rather than
%       sum
% and
%   charge is the charge of the ion (entered by the user);
%   n is a vector of length(isotopeinfo) containing the estimated
%     multiplicity of each element in isotopeinfo.  Only elements with
%     multiple isotopes will have estimates.
%   ncov is the covariance matrix for n, useful for estimating the
%     uncertainties in n (this contains just the non-NaN elements of n)
%   Iout is the output intensity "sparse structure matrix" with the
%     well-fit peaks blanked out.
%
% See also: MSPROF_ISOTOPE_PROFILE.

% Copyright 2009 by Timothy E. Holy

  %% Parse arguments & initialize
  options = struct;
  sProf = [];
  if ~isempty(varargin)
    if isnumeric(varargin{1})
      sProf = varargin{1};
    end
    if isstruct(varargin{end})
      options = varargin{end};
    end
  end
  options = default(options,'peakwidth_mzI',5,'blankwidth_mzI',100,'blank_fit_peaks',nargout > 3 & ~isempty(sProf),...
    'blankfactor',10,'blankmax',1e-3,'zmax',3,'dmzi_max',3,'fit_C2',true,'plot',false);
  % Carbon will be treated specially, prepare for this
  Cindex = strmatch('C',{isotopeinfo.symbol},'exact');
  dC = diff(isotopeinfo(Cindex).isotope_masses(1:2)); % mass increment for charge of 1
  chargeI = zeros(1,options.zmax);

  n_mz = length(I);
  sNorm = 1;
  if ~isempty(sProf)
    sNorm = sum(sProf.^2);
  end

  % Prepare dummy output in case we have to exit early
  charge = nan;
  n_elements = length(isotopeinfo);
  n = nan(1,n_elements);
  ncov = nan(n_elements,n_elements);
  Iout = I;
  Isum_changed = [];
  index_changed = [];

  %% Snip out the base peak
  % We do this early in case we fail to find any isotopologues, because we
  % have to be able to blank this peak to prevent it from coming up again
  width = max(options.blankwidth_mzI,options.peakwidth_mzI);
  indexBase = round(mz2i(mzBase));
  rng0 = (-width:width)+indexBase;
  % Check to make sure we won't "crash" over the edge of the mzIndex
  rng0 = rng0(rng0 > 0);
  rng0 = rng0(rng0 <= n_mz);
  [i1,i2,v] = ssparse_find(I(rng0));
  vtmp = v;
  if ~isempty(sProf)
    vtmp = vtmp.*sProf(i2);
  end
  Ibase = accumarray(i1(:),vtmp(:));
  Ibase(end+1:length(rng0)) = 0;  % pad to full size if there were not enough nonzero elements
  
  %% Blank the base peak
  if options.blank_fit_peaks
    index_changed_c = {rng0};
    Isum_changed_c = cell(1,1);
    Icalc = Ibase(i1)'.*sProf(i2)/sNorm;
    [Iout,Isum_changed_c{1}] = blank_region(Iout,rng0,i1,v,Icalc,options);
  end
  
  % In case of an early exit, prepare outputs
  index_changed = cat(2,index_changed_c{:});
  Isum_changed = cat(1,Isum_changed_c{:});

  %% Determine the charge
  % mzi values for all candidate mz values (+a single 13C for all candidate
  % charge states)
  zvec = (1:options.zmax)/options.zmax;
  indexp1 = round(mz2i(mzBase + zvec*dC));
  if (indexp1(end) > n_mz)
    return
  end
  % Find the candidate charge with highest intensity
  for k = 1:options.zmax
    thisI = I(indexp1(k));
    if ~isempty(sProf)
      chargeI(k) = sum(sProf(thisI.coli) .* thisI.value);
    else
      chargeI(k) = sum(thisI.value);
    end
  end
  if (all(chargeI == 0))
    % Can't determine the charge state, quit
    return
  end
  [mxI,maxIndex] = max(chargeI);
  charge = 1/zvec(maxIndex);
  
  %% Determine which candidate elements have sufficiently-abundant isotopes
  dmzc = cell(1,n_elements);
  ac = cell(1,n_elements);
  for i = 1:n_elements
    if (length(isotopeinfo(i).isotope_masses) > 1)
      dmzc{i} = isotopeinfo(i).isotope_masses(2:end) - ...
	       isotopeinfo(i).base_mass;
      ac{i} = isotopeinfo(i).isotope_abundance(2:end) / ...
	     isotopeinfo(i).isotope_abundance(1);
    end
  end
  l = cellfun(@length,dmzc);
  isotopeIndex = make_counting_index(l);
  dmz = cat(2,dmzc{:})/charge;
  a = cat(2,ac{:});
  % Manually append the 13C2 peak, but initialize it to zero amplitude
  % (we'll fix this later)
  if options.fit_C2
    dmz = [dmz 2*dC/charge];
    a = [a 0];
    isotopeIndex = [isotopeIndex Cindex];
  end
  
  %% Collect the delta m/z into independent regions and snip out the intensity
  % Check to make sure we won't "crash" over the edge of the mzIndex
  % (if that is the case, to simplify & speed the code we just give up)
  if (indexBase - width < 0 || mz2i(mzBase+max(dmz))+width > n_mz)
    return
  end
  % We snip out over a (wide) region of width blankwidth_mzI, but we also supply a
  % mask used for fitting n to restrict the region to peakwidth_mzI within
  % each predicted peak
  mzList = mzBase + dmz;
  mzIndex = mz2i(mzList);
  [regions,regionI] = split_into_contiguous_regions(round(mzIndex),width);
  n_regions = size(regions,2);
  fitmask = cell(1,n_regions);
  for i = 1:n_regions
    rng = regions(1,i):regions(2,i);
    md = mindist(rng,round(mzIndex));
    fitmask{i} = (md <= options.peakwidth_mzI^2)';
  end
  dmzi = zeros(1,n_regions);  % for "registration"
  % Snip out the regions
  Iplus = cell(1,n_regions);
  for i = 1:n_regions
    rng = regions(1,i):regions(2,i);
    [i1,i2,v] = ssparse_find(I(rng));  % deliberately use the original (not that it should matter)
    if ~isempty(sProf)
      v = v .* sProf(i2);
    end
    Iplus{i} = accumarray(i1(:),v(:));
    Iplus{i}(end+1:length(rng)) = 0;  % pad to full size
  end

  %% Use symmetry to increase the likelihood that we're not killing any extra peaks
  % Take the minimimum of Ibase and its flipped version, and
  % use this where there is a significant discrepancy
  Ibasemin = min([Ibase Ibase(end:-1:1)],[],2);
  usemin = Ibase > 2 * Ibasemin;
  Ibase(usemin) = Ibasemin(usemin);

  %% Solve for the number of each element with multiple isotopes
  % We do once for dmzi = 0, then we tweak dmzi, and finally do n again
  % (a cheap version of alternating least squares)
  [n,ncovtmp,Ith] = mspcil_solve_n(Ibase,Iplus,regions,regionI,dmzi,isotopeIndex,mzIndex,a,fitmask);
  if all(isnan(n))
    return;
  end
  if options.fit_C2
    % Manually put in the "quadratic" term for 13C2, leaving out the n
    % coefficient (which is handled by the linear fit). The factor of 2
    % comes from the binomial term (n choose 2).
    a(end) = (isotopeinfo(Cindex).isotope_abundance(2)/isotopeinfo(Cindex).isotope_abundance(1))^2 * max(n(Cindex)-1,0)/2;
    % Generate curves again (this is especially important if there are no
    % other big +2 peaks, e.g., if there are no oxygens or sulfurs)
    [n,ncov,Ith] = mspcil_solve_n(Ibase,Iplus,regions,regionI,dmzi,isotopeIndex,mzIndex,a,fitmask);
  end
  for i = 1:length(Iplus)
    [tmp,dmzi(i)] = register_rigid(Iplus{i},Ith{i},struct('dx_max',options.dmzi_max));
  end
  [n,ncovtmp,Ith] = mspcil_solve_n(Ibase,Iplus,regions,regionI,-dmzi,isotopeIndex,mzIndex,a,fitmask);
  if all(isnan(n))
    return;
  end

  % Fill in NaNs for the undetermined elements
  ncov = nan(n_elements,n_elements);
  ncov(1:length(n),1:length(n)) = ncovtmp;
  n(end+1:n_elements) = nan;
  ncov = ncov / ncov(Cindex,Cindex);  % since we don't have error bars, normalize to carbon
  
  %% If requested, blank out the portions of the spectrum that have been well-described
  if options.blank_fit_peaks
    index_changed_c{n_regions+1} = [];
    Isum_changed_c{n_regions+1} = [];
    for i = 1:n_regions
      rng = regions(1,i):regions(2,i);
      [i1,i2,Imeasured] = ssparse_find(Iout(rng));
      Icalc = Ith{i}(i1)' .* sProf(i2)/sNorm;
      if ~isempty(Imeasured)
        index_changed_c{i+1} = rng;
        [Iout,Isum_changed_c{i+1}] = blank_region(Iout,rng,i1,Imeasured,Icalc,options);
      end
    end
    index_changed = cat(2,index_changed_c{:});
    Isum_changed = cat(1,Isum_changed_c{:});
  end
end

 
%% Utility functions
% Solve for n and produce the "theoretical" (Ithc) intensity profile
function [n,ncov,Ithc] = mspcil_solve_n(Ibase,Iplus,regions,regionsI,dmzi,...
  isotopeIndex,mzIndex,a,fitmask)
  persistent lsqops
  if isempty(lsqops)
    lsqops = optimset(@lsqlin);
    lsqops = optimset(lsqops,'Display','off');
  end
  [Ac,uiI] = msprof_isotope_profile(Ibase,mzIndex,a,isotopeIndex,regions,regionsI,dmzi);
  A = cat(1,Ac{:});
  b = cat(1,Iplus{:});
  fm = cat(1,fitmask{:});
  Asub = A(fm,:);
  if all(Asub(:) == 0)
    n = [];
    ncov = [];
    Ithc = [];
    return
  end
  ns = pinv(Asub)*b(fm); % use pinv to avoid problems with singular matrices
  ns(ns < 0) = 0;
  %ns = lsqnonneg(A(fm,:),b(fm),ns);
  ns = lsqlin(Asub,b(fm),[],[],[],[],zeros(size(A,2),1),[],ns,lsqops);
  n = nan(1,max(uiI));
  n(uiI) = ns;
  if (nargout > 1)
    [~,S,V] = svd(A,'econ');
    Si = diag(1./diag(S));
    ncovs = V*Si*Si*V';
    ncov = nan(length(n),length(n));
    ncov(uiI,uiI) = ncovs;
  end
  if (nargout > 2)
    ns(ns < 0) = 0;  % Don't allow negative numbers of atoms in generating the theoretical intensity profile
    Ithc = cell(size(Iplus));
    for i = 1:length(Iplus)
      Ithc{i} = Ac{i}*ns;
    end
  end
end

% Blank a region
function [Iout,Isum] = blank_region(I,rng,i1,Imeasured,Icalc,options)
% This is an important & tricky function, so a bit of documentation:
%  I: an ssparse matrix, in row format (rows are m/z index)
%  rng: a set of row #s of I to modify
%  i1, Imeasured: these are the values extracted from I, row# (of the
%    subset indexed by rng) and value
%    (note the column number is missing, but we can reconstruct this
%    from I(i1).coli)
%  Icalc: the calculated values at the same points
% On output:
%  Iout: the modified I, blanking the "reasonably well-explained" values
%  Isum: a bit of a misnomer, because the output depends on whether
%    options.usemax is true or false. If false, then Isum is the 
  Iout = I;
  Isum = zeros(length(rng),1);
  % Blank those points that are reasonably well-explained by Icalc
  Imeasured(Imeasured <= options.blankfactor*(Icalc+options.blankmax*max(Icalc))) = 0;
  % Blank those points for which Icalc exceeds a threshold
%   Imeasured(Icalc > options.blank_thresh) = 0;
  breaks_mz = [0 find(diff(i1)>0) length(i1)];
  for indexkill = 1:length(breaks_mz)-1
    % Split out by row #
    thisrng = breaks_mz(indexkill)+1:breaks_mz(indexkill+1);
    % Get the intensities that haven't been set to 0
    thisI = Imeasured(thisrng);
    keepFlag = thisI > 0;
    thisI = thisI(keepFlag);
    i11 = i1(thisrng(1));  % get the relative row #
    indexI = rng(i11);     % get the absolute row #
    if length(Iout(indexI).coli) ~= length(keepFlag)
      warning('Something weird happened')
    end
    % Split the next into two lines for profiling
    tmp = Iout(indexI).coli(keepFlag);
    Iout(indexI).coli = tmp;
    Iout(indexI).value = thisI;
    if options.usemax
      if ~isempty(thisI)
        Isum(i11) = max(thisI);
      else
        Isum(i11) = 0;
      end
    else
      Isum(i11) = sum(thisI);
    end
  end
end
