function stim = default_ud_francesco(type)
% STIMDEFAULT_UD_FRANCESCO: fill in defaults for urine-derived stimuli
  
  stim.category = 'urine-derived';
  stim.identity = '';
  stim.conc_unit = 'relative';
  stim.concentration = NaN;
  stim.strain = 'Balb/c';
  stim.sex = 'F';
  stim.duration = 10;   % 10 s stimulus pulse
  
  switch type
   case 'raw'
    stim.procedure = {};
   case 'filtered'
    stim.procedure = {'filtered'};
   case 'bligh'
    stim.procedure = {'filtered','meth/chlor meth'};
   case 'ODS'
    stim.procedure = {'filtered','meth/chlor meth','C18 elution'};
   case 'WAX'
    stim.procedure = {'filtered','meth/chlor meth','C18 elution','WAX elution'};

  end
  