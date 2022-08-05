function msprof_plot_isotopologues(formula,charge,isotopeinfo,mz,Ismzi,options)
% MSPROF_PLOT_ISOTOPOLOGUES: compare measured and expected abundances of isotopologues
% Syntax:
%   msprof_plot_isotopologues(formula,charge,isotopeinfo,mz,Ismzi,options)
% where
%   formula is a formula vector, corresponding to entries in isotopeinfo,
%     specifying the molecular formula of the candidate ion
%   charge is a positive integer, giving the charge of the ion
%   isotopeinfo contains information provided by "isotopes"
%   mz is a vector of m/z values
%   Ismzi is a matrix of intensities, in the form I(scan#,mzindex)
%   options may have the following fields:
%     mzi_shift (default 0): the shift, in mzindex units, needed to
%       register the measured m/z to the expected m/z for this compound
%     mzi_width (default 10): the number of mzindex values on either side
%       of the expected peak to show
%     plotRawTemporalData (default true): if true, also generates an image
%       plot of Ismzi, which can be useful to insure that one doesn't have
%       two different compounds clouding the analysis
%     figPlot: if specified, uses this figure # for plotting the
%       intensity-vs-mzi data
%     figImage: if specified, uses this figure # for displaying the
%       intensity vs mzi and scan # data (i.e., the temporal plot),
%       assuming plotRawTemporalData is true.
%
% See also: ISOTOPES.

% Copyright 2009 by Timothy E. Holy

  if (nargin < 5)
    options = struct;
  end
  options = default(options,'mzi_shift',0,'mzi_width',20,'plotRawTemporalData',true,'figPlot',-1,'figImage',-1);
  n_mz = length(mz);
  
  [m,a,fstr] = isotopologues(formula,isotopeinfo);
  [m,sortOrder] = sort(m);
  a = a(sortOrder);
  fstr = fstr(sortOrder);
  mzi = interp1(double(mz),1:n_mz,m/charge);  % m/z index associated with each
  rng0 = -options.mzi_width:options.mzi_width;
  % Do the "base" isotopologue and get the m/z profile
  [maxa,maxaIndex] = max(a);
  rng = round(mzi(maxaIndex)+options.mzi_shift) + rng0;
  rng = rng(rng >= 1 & rng <= n_mz);
  ImziBase = full(sum(Ismzi(:,rng),1));
  % Compile all snippets, and compute the expected m/z profile from this
  % isotopologue
  n_isotopologues = length(m);
  mzi_snip = cell(1,n_isotopologues);
  Imzi_th = cell(1,n_isotopologues);
  lbl = ones(1,n_isotopologues);
  %mzi_snip{1} = rng;
  n_snippets = 1;
  a_norm = a / maxa;  % abundance relative to the most-abudant isotopologue
  for i = 1:n_isotopologues
    rng = round(mzi(i)+options.mzi_shift) + rng0;
    rng = rng(rng >= 1 & rng <= n_mz);
    if ~isempty(intersect(mzi_snip{n_snippets},rng))
      % This range overlaps with the previous, just combine into one
      mzi_all = unique([mzi_snip{n_snippets} rng]);
      l = length(mzi_all);
      Ith_all = zeros(1,l);
      thisindx = findainb(mzi_snip{n_snippets},mzi_all);
      Ith_all(thisindx) = Imzi_th{n_snippets};
      thisindx = findainb(rng,mzi_all);
      Ith_all(thisindx) = Ith_all(thisindx) + a_norm(i)*ImziBase;
      mzi_snip{n_snippets} = mzi_all;
      Imzi_th{n_snippets} = Ith_all;
    else
      n_snippets = n_snippets+1;
      mzi_snip{n_snippets} = rng;
      Imzi_th{n_snippets} = a_norm(i)*ImziBase;
    end
    lbl(i) = n_snippets;
  end
  mzi_snip = mzi_snip(1:n_snippets);
  Imzi_th = Imzi_th(1:n_snippets);
  Isnip = cell(1,n_snippets);
  Imzi = cell(1,n_snippets);
  for i = 1:n_snippets
    Isnip{i} = Ismzi(:,mzi_snip{i});
    Imzi{i} = full(sum(Isnip{i},1));
  end
  if ~ishandle(options.figPlot)
    options.figPlot = figure;
  end
  figure(options.figPlot);
  clf
  splitx = SplitAxesEvenly(n_snippets,0.2);
  hax = SplitHoriz(splitx); delete(hax(2:2:end)); hax = hax(1:2:end);
  for i = 1:n_snippets
    axes(hax(i))
    plot(mz(round(mzi_snip{i}-options.mzi_shift)),[Imzi{i}(:) Imzi_th{i}(:)])
    axis tight
  end
  for i = 1:n_isotopologues
    [xn,yn] = data2norm(m(i)/charge,a(i)*max(ImziBase),hax(lbl(i)));
    h = annotation('textarrow',[1 1]*xn,yn*[1.05 1]);
    set(h,'String',fstr{i});
  end
  if options.plotRawTemporalData
    if ~ishandle(options.figImage)
      options.figImage = figure;
    end
    figure(options.figImage);
    clf
    splitx = SplitAxesEvenly(n_snippets,0.2);
    hax = SplitHoriz(splitx); delete(hax(2:2:end)); hax = hax(1:2:end);
    for i = 1:n_snippets
      axes(hax(i))
      imagesc(mz(mzi_snip{i}),1:size(Ismzi,1),Isnip{i});
    end
    set(hax,'YDir','normal')
    ylabel(hax(1),'Scan #')
  end
