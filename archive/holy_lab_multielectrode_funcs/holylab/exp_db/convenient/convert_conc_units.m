function c = convert_conc_units(stimuli,unit)
% CONVERT_CONC_UNITS: change all concentration measurements to a common unit
% Syntax:
%   c = convert_conc_units(stimuli)
%   c = convert_conc_units(stimuli,unit)
% where
%   stimuli is a structure array of stimuli (or a cell array of stimulus
%     structures), which have to have two fields: concentration and
%     conc_unit. conc_unit has to be one of the types accepted by the
%     exp_db, i.e.
%         'none','relative','M','mM','\muM','nM','pM'
%   unit is a string containing the name of the unit that you want to
%     convert to (chosen from one of the above).  The default value is
%     'M'.
% and
%   c is a vector, one element for each of the stimuli, giving the
%     concentration in the chosen unit.

% Copyright 2007 by Timothy E. Holy
  
  if (nargin < 2)
    unit = 'M';
  end

  conc_units = {'none','relative','M','mM','\muM','nM','pM'};
  rel_size =    [NaN      NaN      1  1e-3  1e-6  1e-9 1e-12];
  n_units = length(conc_units);
  conc_unit_to_index = ccu_match_unit(unit,conc_units);
  conc_map = nan(1,n_units);
  if (conc_unit_to_index == 2)
    conc_map(2) = 1;
  elseif (conc_unit_to_index > 2)
    conc_map = rel_size/rel_size(conc_unit_to_index);
  end
  
  n_stimuli = length(stimuli);
  c = nan(1,n_stimuli);
    
  for i = 1:n_stimuli
    if iscell(stimuli)
      this_stimulus = stimuli{i};
    else
      this_stimulus = stimuli(i);
    end
    c(i) = this_stimulus.concentration * ...
	   conc_map(ccu_match_unit(this_stimulus.conc_unit,conc_units));
  end

function conc_unit_index = ccu_match_unit(unit,conc_units)
  conc_unit_index = strmatch(unit,conc_units,'exact');
  if isempty(conc_unit_index)
    fprintf('Valid units:');
    conc_units
    error(['The unit ' unit ' does not match any of these']);
  end
  