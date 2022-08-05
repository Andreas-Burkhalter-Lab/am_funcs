% Example analysis script
% Note this analyzes the different files as if they are separate
% experiments. See below for an alternative if file2 simply contains more
% repeats of the same thing as in file1
bnames = {'file1','file2'};
for i = 1:length(bnames)
  basename = bnames{i};
  fi = imfile([basename '.dat'],[basename '.txt']);
  ip = imphysload(fi);
  % Calculate "Delta f" changes upon stimulation
  df = imstimcalcdf(ip);
  h = fspecial('gaussian',25,8);  % smoothing filter
  smoothfilter = vimage('image',h); % save the definition of the filter
                                    % (this also saves memory)
  dfsmooth = df;   % copy other fields, delete leftovers from previous iter
  dfthumb = df;
  for i = 1:length(df)
    dfsmooth(i).image = imfilter(df(i).image,smoothfilter);
    dfthumb(i).image = imresize(dfsmooth(i).image,0.25,'bilinear',0);
  end
  % Calculate df/f movie
  movie = imphysdf(ip,{[0 1 -1],[-10 -9 -8]});
  % Note the offset starts with 0, so that the "ancestor" of each frame
  % is not offset temporally.
  % Now filter it spatially, and thumbnail it
  for i = 1:length(movie)
    movie(i) = imfilter(movie(i).image,smoothfilter);
    movie(i) = imresize(movie(i),0.25,'bilinear',0);
  end
  % Save results to disk (here's the slow part, because it evaluates the
  % vimages) 
  vimbufresize(50)   % To make calculations more efficient, if the memory
                     % is available
  moviewrite(df,[basename '_df.dat'],dfthumb,[basename '_dfthumb.dat']);
  moviewrite(movie,[basename '_moviethumb.dat']);
  vimagesave([basename '_analyze'],{'raw',ip,'df',df,'dfthumb',dfthumb, ...
                      'movie',movie});
  vimageclear  % You have to do this to clear the globals from memory
end

% If file2 is just basically "more of the same" (more repeats), then you
% can unify the two files by organizing the analysis in the following way:
% for i = 1:nfiles
%    ip{i} = imphysload(...)
%    df{i} = ...
%    movie{i} = ...
% end
% % Now we concatenate
% % Concatenating after analysis means that temporal processing won't
% % accidently grab a frame from the previous file
% ip = [ip{:}];
% df = [df{:}];
% movie = [movie{:}];
% moviewrite(....)
% vimagesave(....)

% To get started on the next stage of analysis at some later point, you'd
% just call "vimageload(filename)" and you're all set to call
% imstimviewdf.
% Here's an example call:
%  imstimviewdf('dfdraw',df,'dfthumb',dfthumb,'movie',movie,'raw',raw,struct('clim',[-.02,.02]))
% This sets the contrast so that 2% changes in DeltaF/F are saturating.
% Try right-clicking on DF/F plots for individual trials.
%
% Note you can also play movies directly:
%  mplayimphys(movie,struct('fps',30,'busymode','drop','showstimulus',1,'clim',[-.02 .02])
