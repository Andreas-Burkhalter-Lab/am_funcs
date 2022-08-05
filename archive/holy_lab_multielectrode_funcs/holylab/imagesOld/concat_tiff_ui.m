inputfiles=UIGetFiles('*.tif', ...
                    'please select tiff files', ...
                    pwd ...
                    );
[filename, pathname] = uiputfile('*.tif', 'Save as ...');
if isequal(filename,0) | isequal(pathname,0)
   disp('User pressed cancel')
else
   fprintf('please wait ... ');
   
   output= fullfile(pathname, filename);
   concat_tiff(inputfiles, output);   
   
   disp('done');
end
                    
