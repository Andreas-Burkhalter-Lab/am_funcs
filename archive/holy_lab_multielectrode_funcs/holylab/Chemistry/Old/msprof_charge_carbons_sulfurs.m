function [charge,nC,nS] = msprof_charge_carbons_sulfurs(mzBase,isotopeinfo,mz2i,Ismzi,varargin)
% MSPROF_CHARGE_CARBONS_SULFURS: determine charge, # of carbons, and # of sulfurs from isotopologues  
%
% Syntax:
%   [charge,nC,nS] = msprof_charge_carbons_sulfurs(mzBase,isotopeinfo,mz2i,I)
%   [charge,nC,nS] = msprof_charge_carbons_sulfurs(mzBase,isotopeinfo,mz2i,I,sProf)
%   [charge,nC,nS] = msprof_charge_carbons_sulfurs(...,options)
% where
%   mzBase is the mean m/z value of the base compound;
%   isotopeinfo is data from ISOTOPES;
%   mz2i is a function that converts m/z values into indices
%   I is the intensity "sparse structure matrix," in row-major ordering
%   sProf (optional) is the profile with respect to scan #;
%   options may have the following fields:
%     peakwidth (default 10): the number of mzi units on either side of the
%       mean m/z to include in the peak
%     zmax (default 4): the maximum charge to consider
%     plot (default false): if true, plots temporal data for isotopologues
% and
%   charge is the charge of the ion (entered by the user);
%   nC is the estimated # of carbons (not necessarily an integer);
%   nS is the estimated # of sulfurs (not necessarily an integer).
%
% The number of carbons and sulfurs is estimated from the abundance of
% +1.0034 and +1.9958 isotopologues. The charge is determined from the
% presence or absence of peaks at +1, +1/2, +1/3, etc.

% Copyright 2009 by Timothy E. Holy

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
  options = default(options,'peakwidth',10,'zmax',4,'plot',false);
  n_mz = length(Ismzi);

  % Gather information about C and S isotopes
  Cindex = strmatch('C',{isotopeinfo.symbol},'exact');
  Sindex = strmatch('S',{isotopeinfo.symbol},'exact');

  % mzi values for all candidate mz values (+a single 13C for all candidate
  % charge states)
  zvec = (0:options.zmax)/options.zmax;
  dC = diff(isotopeinfo(Cindex).isotope_masses(1:2)); % mass increment for charge of 1
  indexp1 = round(mz2i(mzBase + zvec*dC));
  indexBase = indexp1(1);
  % Find the candidate charge with highest intensity
  chargeI = zeros(1,options.zmax);
  for k = 1:options.zmax
    thisI = Ismzi(indexp1(1+k));
    if ~isempty(sProf)
      chargeI(k) = sum(sProf(thisI.coli) .* thisI.value);
    else
      chargeI(k) = sum(thisI.value);
    end
  end
  [mxI,maxIndex] = max(chargeI);
  charge = 1/zvec(maxIndex+1);
  
  % Determine m/z of 13C and 34S isotopologues
  dC = diff(isotopeinfo(Cindex).isotope_masses(1:2))/charge;
  dS = diff(isotopeinfo(Sindex).isotope_masses(1:2))/charge;
  mzi = round(mz2i([dC 2*dC dS]+mzBase));
  % Check to make sure we won't "crash" over the edge of the mzIndex
  % (if that is the case, to simplify & speed the code we just give up)
  if (indexBase - options.peakwidth < 0 || max(mzi)+options.peakwidth > n_mz)
    nC = nan;
    nS = nan;
    return
  end
  
  rng0 = -options.peakwidth:options.peakwidth;
  % Calculate the total intensity of the base peak
  [i1,i2,v] = ssparse_find(Ismzi(rng0+indexBase));
  if ~isempty(sProf)
    basepeakI = sum(v .* sProf(i2));
  else
    basepeakI = sum(v);
  end
  % Estimate the number of carbons from the +1 peak
  [i1,i2,v] = ssparse_find(Ismzi(rng0+mzi(1)));
  if ~isempty(sProf)
    Ip1 = sum(v .* sProf(i2));
  else
    Ip1 = sum(v);
  end
  nC = (Ip1/basepeakI) / (isotopeinfo(Cindex).isotope_abundance(2)/isotopeinfo(Cindex).isotope_abundance(1));
  % Estimate the number of sulfurs from the +1.9958 peak
  rng1 = mzi(3)-options.peakwidth:min(mzi(3)+options.peakwidth,round(sum(mzi(2:3))/2)); % be sure to stop before 13C2 peak
  [i1,i2,v] = ssparse_find(Ismzi(rng1));
  if ~isempty(sProf)
    IS = sum(v .* sProf(i2));
  else
    IS = sum(v);
  end
  nS = (IS/basepeakI) / (isotopeinfo(Sindex).isotope_abundance(2)/isotopeinfo(Sindex).isotope_abundance(1));
  % Plot data for user
  if options.plot
    regions = split_into_contiguous_regions([indexBase mzi],2*options.peakwidth);
    figure
    n_regions = length(regions);
    xsplit = SplitAxesEvenly(n_regions,0.2);
    hax = SplitHoriz(xsplit); delete(hax(2:2:end)); hax = hax(1:2:end);
    for i = 1:n_regions
      rng = regions(1,i):regions(2,i);
      axes(hax(i))
      imagesc(mz(rng),1:size(Ismzi,1),Ismzi(:,rng));
      set(gca,'YDir','normal')
    end
    xlabel(hax(1),'Scan #')
  end