function ui_imstimplay(timerange)
   if(isunix)
      file_filter='*';
   else
      file_filter='*.*';
   end
   
   [filename, pathname] = uigetfile(file_filter, 'Pick a data file');
   if(filename==0)
      disp('user cancelled');
      return;
   end

   data_filename  =fullfile(pathname, filename);
   
   [filename, pathname] = uigetfile(file_filter, 'Pick a description file');
   if(filename==0)
      disp('user cancelled');
      return;
   else
      header_filename=fullfile(pathname, filename);
   end
   
   %  data_filename  ='/home/jason/raw_test/raw_test';
   %  header_filename='/home/jason/raw_test/raw_test.txt';

   if(nargin==0)
      timerange=[-12 80];
   end
   
   ip=imphysfrom2d(data_filename, header_filename);
   imstimplay(ip, struct('trange', timerange ) );
