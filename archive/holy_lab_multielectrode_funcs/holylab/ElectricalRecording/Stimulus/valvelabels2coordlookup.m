function [coordlookup,uid] = valvelabels2coordlookup(valvelabels)
% VALVELABELS2COORDLOOKUP: generate "coordinates" corresponding to stimulus identity
% Syntax:
%    [coordlookup,identity] = valvelabels2coordlookup(valvelabels)
% where
%   valvelabels is a cell array of strings describing your stimulus
%     contents
% and
%   coordlookup is a matrix of the form described in stim2coord
%   identity is a cell array of strings describing the coordinates (in
%     order)
%
% See also: STIM2COORD.

% Copyright 2008 by Timothy E. Holy

  nvalves = length(valvelabels);
  stimuli = parse_valvelabels(valvelabels);
  conc = nan(1,nvalves);
  id = cell(1,nvalves);
  for i = 1:nvalves
    id{i} = [stimuli{i}.category ' ' stimuli{i}.identity];
    conc(i) = convert_conc_units(stimuli{i});
    if strcmp(stimuli{i}.category,'urine-derived')
      conc(i) = stimuli{i}.concentration;
    end
  end
  
%   issl = false(1,nvalves);
%   for i = 1:nvalves
%     [issl(i),cIndex] = issteraloids(valvelabels{i});
%     if issl(i)
%       id{i} = valvelabels{i}(1:cIndex-1);
%       conc(i) = concstring2molar(valvelabels{i}(cIndex:end));
%     elseif any(cellfun(@(s) ~isempty(findstr(s,lower(valvelabels{i}))),negctrls))
%       id{i} = 'neg';
%     else
%       id{i} = valvelabels{i};
%     end
%   end

  [uid,tmp,idnum] = unique(id);

  nuid = length(uid);
  coordlookup = zeros(nuid,nvalves+1);  % +1 is for flush
  clabel = agglabel(idnum);
  for i = 1:nuid
    this_conc = conc(clabel{i});
    if any(isnan(this_conc))
      % Just set this coordinate to 1, because we don't have concentration
      % info
      coordlookup(i,clabel{i}+1) = 1;
    else
      % Arrange the concentrations in order
      [uc,tmp,uorder] = unique(this_conc);
      coordlookup(i,clabel{i}+1) = uorder/max(uorder);
    end
  end
end
