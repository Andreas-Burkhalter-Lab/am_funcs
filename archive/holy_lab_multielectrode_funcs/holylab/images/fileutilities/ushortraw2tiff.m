function ushortraw2tiff(filename, options)

   d=dir(filename);
   filesize=d.bytes;
   
   if(nargin==1)
      imageW=1360;
      imageH=1024;
      reverse_row_order=0;
   else
      imageW=options.imageW;
      imageH=options.imageH;
      reverse_row_order=options.reverse_row_order;
   end
   
   if(filesize<imageW*imageH*2)
      error('raw file is incorrect');
   end
   
   fd=openlfs(filename);
   if(fd<0)
      error('failed to open file');
   end
   
   tName=[' %d %% done'];
   figProgress=waitbar(0,['converting ' filename  ', please wait...'], 'name', sprintf(tName, 0));
   
   tiff_filename=[filename '.tif'];
   nRowGot=0;
   frame=readint16lfs(fd, imageW, nRowGot+[0, imageH-1], 0);
   frame=frame';
   if(reverse_row_order) frame=frame(end:-1:1,:); end
   nRowGot=nRowGot+imageH;
   imwrite(uint16(frame), ...
        tiff_filename, ...
           'tif', ...
           'writemode', 'overwrite', ...
           'compression','none',...
           'description', 'first frame');
   
   nFramesRemain=floor((filesize-imageW*imageH*2)/(imageW*imageH*2));
   for idx=1:nFramesRemain
      frame=readint16lfs(fd, imageW, nRowGot+[0, imageH-1], 0);
      frame=frame';
      if(reverse_row_order) frame=frame(end:-1:1,:); end
      nRowGot=nRowGot+imageH;   
      imwrite(uint16(frame), ...
              tiff_filename, ...
              'tif',...
              'writemode', 'append', ...
              'compression', 'none',...
              'description', [num2str(idx+1) 'th frame']);
      
      % show progress:
      if(ishandle(figProgress) & mod(idx, 10) == 0) 
         set(figProgress, 'name', sprintf(tName, round(idx/nFramesRemain*100))); 
         waitbar(idx/nFramesRemain, figProgress);
         % drawnow;
      end
   end
   
   closelfs(fd);

   if(ishandle(figProgress)) close(figProgress); end
   
