function fix_stim_seq()
% fix the "stimulus sequence" field in a merec file according to related stim.vlv file
% @usage: fix_stim_seq()
% @note: use this function to fix header of merec data created by old merec executable.
% @see: is_new_merec()
% 
   [tfilename, tpath]=uigetfile('*.merec', 'pls choose a .merec file');
   merec_file=[tpath filesep tfilename];
   
   vlv_file=[tpath filesep 'stim.vlv'];
   if(exist(vlv_file, 'dir') | ~exist(vlv_file, 'file'))
      vlv_file=[tpath filesep 'old analysis' filesep 'stim.vlv'];
      if(exist(vlv_file, 'dir') | ~exist(vlv_file, 'file'))
         [tvlvfilename, tvlvpath]=uigetfile('*.vlv', 'pls choose a .vlv file');
         vlv_file=[tvlvpath filesep tvlvfilename];
      end
   end
   
   [tt1,basename,tt2] = fileparts(merec_file);
   cmd=['grep "' basename ' "  "' vlv_file '" | cut -d \{ -f 2 | cut -d \} -f 1'];
   [status, output]=system(cmd);
   stim_seq=eval(['[' output ']']);
   tt=zeros(2, size(stim_seq,2)/2);
   tt(:)=stim_seq(:);
   stim_seq=tt';
   ttHeader=readheader(merec_file);
   stim_seq(:,2)=stim_seq(:,2)/ttHeader.scanrate; % scan number to second
   
   strStim = num2str(stim_seq(1,:));
   for idx=2:size(stim_seq, 1) 
      strStim=[strStim ';' num2str(stim_seq(idx,:))];
   end
   
   % strStim
   
   cmd=['echo -n "' strStim '" | setheader -k "stimulus sequence" -f ' merec_file];
   [status, output]=system(cmd);
   
   
