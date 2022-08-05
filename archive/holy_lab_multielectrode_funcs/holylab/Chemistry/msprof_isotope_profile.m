function [A,uelementIndex] = msprof_isotope_profile(Ibase,mzIndex,a,elementIndex,regions,regionI,dmzI)
% MSPROF_ISOTOPE_PROFILE: create matrices facilitating isotopologue fitting
%
% This function is valid only for single-substitution isotopologues.  For
% small molecules, these are the most abundant, and thus provide a basis for
% fast fitting because one can use linear algebra.  For more complete
% distributions, see ISOTOPOLOGUES.
%
% Syntax:
%   [A,uelementIndex] =
%      msprof_isotope_profile(Ibase,mzIndex,a,elementIndex,regions,regionI,dmzI)
% where
%   Ibase is the m/z profile of the "base" ion; must be sampled
%     symmetrically around the peak (i.e., from -peakwidth:peakwidth
%     relative to the peak)
%   mzIndex is a vector of the m/z index associated with each possible
%     single-substitution isotopologue (can be fractional)
%   a is a vector of the abundance (relative to the base) of the
%     substituted isotope (not including the multiplicity in the
%     molecular formula---that's what we're trying to deduce!)
%   elementIndex is a vector providing an integer identifier for the type
%     of atom being substituted (note one can thereby support multiple
%     isotopes for each element)
%   regions, regionsI are in the output format of
%     split_into_contiguous_regions, and describe the contiguous ranges
%     of m/z indices that are going to be fit (often corresponding to the
%     +1 isotopologues, the +2 isotopologues, etc.)
%   dmzI is a vector of length 1-by-n_regions that specifies the amount
%     to shift the base profile by in each region (allows for
%     "registration")
% and
%   A is a cell array of length n_regions, each element with format
%     described below;
%   uelementIndex is the list of unique element indices.
%
% Suppose n is the vector of element abundances, and un =
% n(uelementIndex).  Then in the ith region, one generates the predicted
% m/z profile by
%    Ith = A{i} * un(:);
%
% An example of usage is found in msprof_charge_isotopes_linear.
% See also:  MSPROF_CHARGE_ISOTOPES_LINEAR, ISOTOPOLOGUES.
  
% Copyright 2009 by Timothy E. Holy
  
  [uelementIndex,tmp,elementIndex] = unique(elementIndex);
  n_elements = length(uelementIndex);
  n_peaks = length(mzIndex);
  len = diff(regions,1,1)+1;
  halfwidth = (length(Ibase)-1)/2;
  n_regions = size(regions,2);
  A = cell(n_regions,1);
  for i = 1:n_regions
    A{i} = zeros(len(i),n_elements);
  end
  for i = 1:n_peaks
    thisRegion = regionI(i);
    thisElement = elementIndex(i);
    t = mzIndex(i) + dmzI(thisRegion) - regions(1,thisRegion) - halfwidth + 1;
    v = template_copies(Ibase,t,a(i),[1 len(thisRegion)]);
    A{thisRegion}(:,thisElement) = A{thisRegion}(:,thisElement) + v;
  end
end
