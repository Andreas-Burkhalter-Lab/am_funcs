function raw2tiff(options)

default_options('width', 1004);
default_options('height', 1002);

rawfiles=UIGetFiles('*', ...
                    'please select raw data files' ...
                    );

for idx=1:length(rawfiles)
   ushortraw2tiff(rawfiles{idx}, struct('imageW', options.width, ...
      'imageH', options.height, 'reverse_row_order', 1));
end
