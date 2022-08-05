function result=gen_xdb_entry_behavior(dirs, options)
% this func generates an entry for a behavior experiment phase, of which
% files were saved under dirs
% SYNTAX: result=gen_xdb_entry_behavior(dirs, options)
% PRE:
%    dirs: a cell array of strings

   default_options('investigator', {getenv('USER')});

   result.version=1;
   
   result.investigator=options.investigator;
   
   if(~iscell(dirs))
       dirs={dirs};
   end
   
   dates=datenum(dirs, 'yyyymmdd');
   dates=sort(dates);
   
   result.start_date=datevec(dates(1));
   result.start_date=result.start_date(1:3);
   
   result.tags={'singing'};
   
   % now stimuli field
   if(isequal(result.investigator, {'jason'}))
      stimuli={};
      
      stim.category='urine_whole';
      stim.identity='M';
      stim.concentration=0.3162;
      stim.conc_unit='fold';
      stim.duration=180; % in sec
      stim.usr_tag='';
      stim.strain='BALB/cJ';
      stim.sex='M';
      
      stim4mixMale=[stim stim stim];
      concentrations=num2cell(10.^[-0.5 -1 -1.5]);
      [stim4mixMale.concentration]=deal(concentrations{:});
      
      stimuli{end+1}=stim; % pure male urine
      
      stim.identity='F';
      stim.sex='F';
      stimuli{end+1}=stim; % pure female urine
      
      stim4mixFemale=[stim stim stim];
      [stim4mixFemale.concentration]=deal(concentrations{end:-1:1});
      
      identity=1;
      for i=1:length(stim4mixFemale)
         for j=1:length(stim4mixMale) 
            stimuli{end+1}=[stim4mixFemale(i) stim4mixMale(j)]; % mixture
            stimuli{end}(1).identity=num2str(identity);
            stimuli{end}(2).identity=num2str(identity);
            identity=identity+1;
         end
      end
      
      % now blank stim
      stim.substance='';
      stim.identity='0';
      stim.concentration=0;
      stim.conc_unit='fold';
      stim.duration=180; % in sec
      stim.usr_tag='blank';
      stim.strain='';
      stim.sex='';
      
      stimuli{end+1}=stim;
      
      switch(options.phase)
          case 'a', result.stimulus=stimuli([12 3:11 2 1]); % in the order of identity: 0-9, f,m
          case 'b', result.stimulus=stimuli([12 2 1]); 
          case 'd', result.stimulus=stimuli([12 3:11 2 1]);
          case 'e', result.stimulus=stimuli([12 2 1]);
          case 'f', result.stimulus=stimuli([12 3:11 2 1]);
          otherwise, error('unknow phase');
      end
   else
      error('please update the code to reflect your protocol'); 
   end
   
   result.data_locations=strcat('$EXP_DATA_ROOT$/', dirs);
   
   % comment field:
   if(isequal(result.investigator, {'jason'}))
      switch(options.phase)
          case {'a', 'b', 'd', 'e', 'f'} 
             result.comment=['phase ' options.phase]; 
          otherwise, error('unknow phase');
      end
   else
      error('please update the code to reflect your protocol');
   end
   
   
