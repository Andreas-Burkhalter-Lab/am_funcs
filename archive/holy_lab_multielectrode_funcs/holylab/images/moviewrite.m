function moviewrite(varargin)
% MOVIEWRITE: write a movie to disk
%
% This function allows you to write one or more movie files.
%
% Syntax:
%   moviewrite(ip,filename)
%   moviewrite(ip1,filename1,ip2,filename2) (see "Efficiency note" below)
%   moviewrite(...,'fake')
% where
%   ip are IMPHYS structure arrays;
%   filename are strings containing the name of the file.
%
% Note that there are no particular constraints on the movie file
% format, i.e., you can combine frames of different sizes, etc, into one
% file. Note that concatenation:
%
%   moviewrite([ip1 ip2],filename)
% 
% is a perfectly reasonable way to get two types/versions of a movie into
% the same file.
%
% If the supplied movies are vimages, then moviewrite pushes the "imload"
% method onto the vimage stack so that the movie can be easily loaded in
% the future. See VIMAGESAVE to save the vimage data.
%
% If the final argument is the string 'fake', then the images are not
% actually written to disk---or for vimages, evaluated.  However, the
% file position indicators are incremented appropriately as if writing
% has occurred.  This is useful if you need to "re-run" an analysis to
% change the data saved by vimagesave, as long as the movie itself does
% not need re-writing.
%
% Efficiency note: suppose you have two versions of a movie, "smooth" and
% "thumb" (a "thumbnail" movie) and "thumb" is based on "smooth," in the
% sense that "thumb" is specified by vimages that calculate from
% "smooth." You can take advantage of the fact that the frames of
% "smooth" will be buffered to give you a head start on calculating
% "thumb." In particular, if you call
%     moviewrite(smooth,'smooth.dat',thumb,'thumb.dat')
% then the evaluated "smooth" images serve as a starting point for
% "thumb" images, whereas if you call
%     moviewrite(smooth,'smooth.dat')
%     moviewrite(thumb,'thumb.dat')
% then each frame of "thumb" must first re-calculate "smooth."
%
% If, however, there is no relationship between the movies, then there is
% no reason to specify the two in a single call. Indeed, if temporal
% processing is involved, it may be a disadvantage.
%
% Finally, similar remarks apply to concatenated movies, [ip1 ip2],
% written to a single file. In calling MOVIEWRITE, you may want to
% arrange the frames in an order which optimizes buffer utilization. You
% can of course arrange them in any sequence you wish for playback.
%
% See also: VIMAGESAVE, IMWRITE, VIMAGELOAD.
  
% Copyright 2005 by Timothy E. Holy
  
  % Parse input arguments
  npairs = floor(length(varargin)/2);
  movie = varargin(1:2:2*npairs);
  filename = varargin(2:2:2*npairs);
  mode = '';
  if (length(varargin) > 2*npairs)
    mode = varargin{end};
  end

  % Determine the lengths of the movies
  len = zeros(1,length(movie));
  for i = 1:length(movie)
    len(i) = length(movie{i});
  end
  maxlen = max(len);

  % Open output files
  if strcmp(mode,'fake')
    fid = zeros(1,length(filename));  % This will hold the fake file
                                      % position
    [str,maxsize,machfmt] = computer;
    fh = filehandle;
    for i = 1:length(filename)
      fh(i) = filehandle('filename',filename{i},...
                         'machfmt',lower(machfmt));
    end
  else
    for i = 1:length(filename)
      [fid(i),msg] = fopen(filename{i},'w');
      if (fid(i) < 0)
        error(['moviewrite:fopen: ' msg]);
      end
    end
  end
  
  % Create name for progress report
  filestr = [filename; repmat({'; '},1,length(filename))];
  filestr = [filestr{:}];
  filestr = filestr(1:end-2);  % Strip off last semicolon
  if ~isempty(mode)
    filestr = [filestr ' (' mode ')'];
  end
  
  % Start writing frames; work through all the movies in parallel
  counter = 1;
  while (counter <= maxlen)
    if (~strcmp(mode,'fake') || ~mod(counter-1,ceil(maxlen/10)) || counter==maxlen)
      progress_bar(struct('progress',counter,'max',maxlen,...
        'what',['Saving ' filestr]));
    end
    for i = 1:length(filename)
      if (counter <= length(movie{i}))
        frame = movie{i}(counter);
        if strcmp(mode,'fake')
          [fileinfo,fid(i)] = imwrite(fid(i),frame.image,mode);
          fileinfo.fh = fh(i);
        else
          fileinfo = imwrite(fid(i),frame.image);
        end
        fileinfo = copyfields(frame,...
              {'stacknum','stacktime','stimulus'},fileinfo);
        if isa(frame.image,'vimage')
          push(frame.image,'imload',fileinfo);
        end
      end
    end
    counter = counter+1;
  end

  fprintf('Don''t forget to call vimagesave\n');
  