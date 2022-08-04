%----------------------------------------------------------------------------------------------------------
%-- GetHTBInfo.m: This function extracts various useful pieces of information from the TEMPO HTB headers
%--	and returns this database information in a cell array of strings, info_text.
%--	GCD, 12/31/99
%----------------------------------------------------------------------------------------------------------

% 3/15/16 AM changed capitalization make compatible with new matlab
% versions

function	[return_value, info_text] = GetHTBInfo(PATH,FILE)

	TEMPO_Defs;		%reads in some definitions that we need

	l = length(FILE);
   if (FILE(l-3:l) == '.htb')	% .htb extension already there
		filename = [PATH FILE];   %the HTB data file
      logfile = [PATH FILE(1:l-4) '.log'];   %the TEMPO log file
   else	%no extension in FILE, add extensions
		filename = [PATH FILE '.htb'];   %the HTB data file
      logfile = [PATH FILE '.log'];   %the TEMPO log file
   end

	fid = htbOpen(filename); % Open the HTB file % % 3/15/16 AM changed capitalization make compatible with new matlab versions
   if (fid == -1)		%file could not be opened
      return_value = -1;
      return;
   end
   
   info_text = [];
   
   %read in protocol type and add to cell array
   protocol_info = textread(logfile, '%s', 3);
   protocol_name = protocol_info{3};
   t = sprintf('Protocol Type = %s', protocol_name);
   info_text{length(info_text)+1} = t;
   
   ndbs = htbCount(fid); % Get number of databases in it % % 3/15/16 AM changed capitalization make compatible with new matlab versions
   h = htbGetHd(fid, 1);   % HTB header #1 %%%% 3/15/16 AM changed capitalization
	t = sprintf('File: %s', filename);
   info_text{length(info_text)+1} = t;		%add string to cell array
	t = sprintf('# of Databases = %d', ndbs);
   info_text{length(info_text)+1} = t;		%add string to cell array
   
	%print out some useful information
	for database = 1:ndbs               % Loop though each one...
	   h = htbGetHd(fid, database);   % HTB header %%%% 3/15/16 AM changed capitalization
      
      t = sprintf('Database #%d: %s, # of channels = %d, # of epochs = %d', database, h.title, h.nchannels, h.sweep);
     	info_text{length(info_text)+1} = t;
		t = sprintf('Collected on %s using %s', h.date, h.pro_file);
   	info_text{length(info_text)+1} = t;		%add string to cell array

	   hertz = h.speed_units / h.speed;    % see p366 in TEMPO v9 manual
   	binwidth = (h.skip + 1) / hertz;
	   epoch_start = (-h.offset) * binwidth;
   	epoch_period = h.period * binwidth;
   	t =sprintf('    Sample Rate = %0.3f (Hz),  Time Window = %0.3f -> %0.3f (sec)', ...
      	(1/binwidth), epoch_start, epoch_start + epoch_period);
     	info_text{length(info_text)+1} = t;
	end
   
	err = htbClose(fid);   % Close HTB file  %%%% 3/15/16 AM changed capitalization
   return_value = 1;		%indicates went OK
return;