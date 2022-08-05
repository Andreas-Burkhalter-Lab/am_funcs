function suc=trylock_file(filename, options)
% try to lock a file
% SYNTAX:
%   suc=trylock_file(filename, options)
% PRE:
%   filename: the file to lock
%   options: a struct. Currently 2 fields are supported:
%      .retries=2: # of retries. Put 0 here for no retry.
%      .max_sleep_time=5: after each try, max time to sleep in sec.
% NOTE:
%   the lock is advisory and currently works on only linux.
   
   default_options('retries', 2); % i.e. try 3 times
   default_options('max_sleep_time', 5);
   
   while(true)
      suc=trylock_file_noretry(filename);
      if(suc) return; end
      if(options.retries<=0) return; end
      sec2sleep=options.max_sleep_time*rand;
      sleep_g(sec2sleep);
      options.retries=options.retries-1;
   end

   
function suc=trylock_file_noretry(filename)
   cmd=['dotlockfile -p -r 0 -l ' filename '.lock'];
   [st, tt]=system(cmd);
   suc=st==0;
      