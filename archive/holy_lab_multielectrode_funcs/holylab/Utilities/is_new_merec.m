function result=is_new_merec(merec_file)
% @usage:  result=is_new_merec(merec_file)
% @note: return 1 if the file was created by new merec executable
%        return 2 if unknown(when no stimulus sequence is present in the header),
%        return 0 if the file was created by old merec executable
% 
   ttHeader=readheader(merec_file);
   stim_from_file=parse_stim_seq(ttHeader);
   if(isempty(stim_from_file))
      result=2;
      return;
   end
   
   [tpath,basename,tt] = fileparts(merec_file);
   if(~isempty(tpath))
      tpath=[tpath filesep];
   end
   
   vlv_file=[tpath 'stim.vlv'];
   if(exist(vlv_file, 'dir') | ~exist(vlv_file, 'file'))
      vlv_file=[tpath 'old analysis' filesep 'stim.vlv'];
      if(exist(vlv_file, 'dir') | ~exist(vlv_file, 'file'))
         [tvlvfilename, tvlvpath]=uigetfile('*.vlv', 'pls choose a .vlv file');
         vlv_file=[tvlvpath filesep tvlvfilename];
      end
   end
   
   cmd=['grep "' basename ' "  "' vlv_file '" | cut -d \{ -f 2 | cut -d \} -f 1'];
   [status, output]=system(cmd);
   stim_seq=eval(['[' output ']']);
   tt=zeros(2, size(stim_seq,2)/2);
   tt(:)=stim_seq(:);
   stim_from_vlv=tt';
   stim_from_vlv(:,2)=stim_from_vlv(:,2)/ttHeader.scanrate; % scan number to second
   
   nrow=min(size(stim_from_file, 1), size(stim_from_vlv, 1));
   stim_from_file=stim_from_file(1:nrow, :);
   stim_from_vlv =stim_from_vlv(1:nrow, :);
   
   ttMaxTimeDiff=max(abs(stim_from_vlv(:,2)-stim_from_file(:,2)));
   
   if(ttMaxTimeDiff>1)
      result=0;   
   else
      result=1;
   end
