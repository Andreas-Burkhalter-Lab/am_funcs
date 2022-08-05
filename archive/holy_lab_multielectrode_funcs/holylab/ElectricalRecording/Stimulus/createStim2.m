function createStim2()
% CREATESTIM2: create stimulus files
% syntax:
%    createStim2()
% obsoletes: createStim()
% 
% @history: since 20031001:
%    20031001: solve the problem "error while loading shared libraries: libborqt-6.9.0-qt2.3.so"
% 
   
   if(exist('/usr/local/kylix3/bin', 'dir')==7)
      kylixlibpath='/usr/local/kylix3/bin';
   else
      % [tt, kylixlibpath]=system('echo -n $HOME'); % this work
      kylixlibpath=getenv('HOME'); % a better way
      kylixlibpath=[kylixlibpath filesep 'kylix3/bin'];
   end

   [status, output]=system(['env LD_LIBRARY_PATH=' kylixlibpath ' createStim2 -o --matlab']);

   stim=eval(output);
   if(~isempty(stim))
      figure;
      stairs(stim(:,2),stim(:,1));
      axis tight
      title(['Here is the last stimulus sequence created']);
   end
