function view3d(src, options)
% view or take a snapshot of a stack volume
% USAGE: view3d(imagine, stack, outdir, action, rotation)
% PRE:
%    src: a .imagine file or a stackmm object or a struct
%         (the stuct has 3 field: data, um_per_pixel and imagine_filename. The last one is used only to
%          deduct the path and name of the uvf file)
%    options: a struct w/ possible fields:
%       stack: the stack # (NOTE: it's numbered from 1)
%       outdir: where the tmp file and snapshot go. Use '' for the default.
%       action: one of 'capture', 'view', 'conversion'
%       rotation: a 1x3 vector specifying the degrees to rotate
%       one_d_file: 1d transfer func setting file

if(nargin==1)
   options=struct; 
end
options=default(options, 'stack', 1, 'outdir', '', 'action', 'view', 'rotation', [0 0 0], 'one_d_file', '');

stack=options.stack;
outdir=options.outdir;
action=options.action;
rotation=options.rotation;
one_d_file=options.one_d_file;

hasData=false;
if(isobject(src))
   imagine=src.filename;
   imagine=replace_extension(imagine, '.imagine'); % add the ext if nec
   smm=src;
elseif(isstruct(src))
   imagine=src.imagine_filename;
   hasData=true;
else
   imagine=src;
   smm=stackmm(imagine);
end

% imagine='/usr/lab/home/jason/tmp4_3d/2006-01-31-27umbeads.imagine';
cam=replace_extension(imagine, '.cam');
if(isempty(outdir))
   outdir=replace_filename(imagine, ''); % the parent dir of imagine
end

rawfilename=fullfile(outdir, ...
   [replace_extension(replace_parent_dir(cam, ''), ''), ...
   sprintf('.%04d', stack), ...
   '.cam']);

if(hasData)
   h.eff_dim=size(src.data);
   h.prec='uint16'; % todo: hardcoded
   h.um_per_pixel=src.um_per_pixel;
else
   h=smm.header;
   if(~isfield(h, 'eff_width'))
      h.eff_dim=smm.size;
   end
   if(~isfield(h, 'um_per_pixel'))
      h.um_per_pixel=[h.um_per_pixel_xy h.um_per_pixel_xy diff(h.piezo_start_stop)/(h.frames_per_stack-1)];
   end
end

uvffilename=replace_extension(rawfilename, '.uvf');

if(~fileexist(uvffilename))
   fid=fopen(rawfilename, 'w');
   if(fid<0)
      error('failed to write raw stack data');
   end
   if(hasData)
      data=src.data;
   else
      data=smm(:,:,:, stack);
   end
   % data=permute(data, [2 1 3]);
   count=fwrite(fid, data, h.prec);
   if(count~=numel(data))
      error('failed to write raw stack data: not all data was written');
   end
   fclose(fid);
end

scriptname=replace_extension(uvffilename, '.script');

if(fileexist(uvffilename))
   script=sprintf('open "%s"\n', uvffilename);
else
   script=sprintf('open "%s" "%s" %d,%d,%d,%f,%f,%f\n', rawfilename, uvffilename, ...
      h.eff_dim(1), h.eff_dim(2), h.eff_dim(3), h.um_per_pixel(1), h.um_per_pixel(2), ...
      h.um_per_pixel(3));
end

if(~isequal(action, 'conversion'))
   script=[script ...
      sprintf('rotateX %d\nrotateY %d\nrotateZ %d\n', rotation(1), rotation(2), rotation(3))];
   if(~isempty(one_d_file))
      script=[script ...
         sprintf('open1d "%s"\n', one_d_file)];
   end
end

if(isequal(action, 'capture'))
   script=[script ...
      sprintf('capturesingle "%s"\n', replace_extension(uvffilename, '.png'))];
elseif(isequal(action, 'view'))
   script=[script ...
      sprintf('stayopen\n')];
end

fid=fopen(scriptname, 'w');
if(fid<0)
   error('failed to write script file');
end
fwrite(fid, script, 'char');
fclose(fid);

cmd=sprintf('env LD_LIBRARY_PATH="" /usr/bin/imagevis3d_g -script "%s"', scriptname);
system(cmd);
