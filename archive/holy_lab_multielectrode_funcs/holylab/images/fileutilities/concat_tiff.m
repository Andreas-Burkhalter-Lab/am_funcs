function concat_tiff(inputs, output)
%  concat_tiff: concatenate tiff files into one multi-tiff file
%  usage:
%    concat_tiff(inputs, output)
%  pre: 
%     inputs: a cell array of filenames
%     output: the generated tiff file's name
   isFirstFrame=1;
   for idx=1:length(inputs)
      curFile=inputs{idx};
      iminfs =  imfinfo(curFile);
      for idxFrame=1:length(iminfs)
         data = imread(curFile, idxFrame);
         if(isFirstFrame)
            isFirstFrame=0;
            writemode='overwrite';
         else
            writemode='append';
         end
         
         imwrite(data, ...
                 output, ...
                 'tif', ...
                 'writemode', writemode, ...
                 'compression','none',...
                 'description', [curFile ', frame idx (1-based) ' num2str(idxFrame)]);
         
      end
   end
