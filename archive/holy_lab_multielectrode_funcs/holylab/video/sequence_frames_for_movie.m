function sequence_frames_for_movie(basenamein,sequence,basenameout,options)
% sequence_frames_for_movie: arrange frames in a sequence for conversion via ffmpeg
%
% Given a directory with a series of image files with the same base name,
% this creates a set of symbolic links that put the frames into a
% specified order. These links can be used (e.g., with ffmpeg) to make a
% movie.
%
% Let's start with an example: you have a directory with frames labelled
%     frame1.png
%     frame2.png
%     ...
%     frame40.png
% Suppose you want to make a movie playing at 2 frames per second. You want
% to make a movie that shows frame 1-20 in order, then pauses on frame 20
% for 5 seconds, then shows frames 1-20 again but proceeds immediately
% onward to frames 21-40. You might arrange this in the following way:
%    sequence = [1:20 repmat(20,1,10) 1:40];
% Here, 20 gets repeated 10 times because (5 sec)*(2 frames/sec) = 10
% frames.
% Then you'd make the links this way:
%   sequence_frames_for_movie('frame',sequence,'framelink')
%
% Here's the general syntax:
%   sequence_frames_for_movie(basenamein,sequence,basenameout,options)
% where
%   basename is the base name of your movie frames (note: if your input
%     frames use a zero-padding convention, e.g., frame 3 has filename
%     'frame003.png', then set options.pad_with_zero as described below)
%   sequence is a vector specifying the order in which frames appear
%   basenameout is the base name of the symbolic links that will be
%     created. The file names will be zero-padded automatically.
%   options is a structure which may have the following fields:
%     pad_with_zero (default 0): if greater than zero, the input file names
%       are zero-padded up to pad_with_zero digits.
%     overwrite_automatically (default false): if true, any files that
%       begin with basenameout are first deleted. Use with caution!
%     extension: a string like '.png' that specifies the extension of your
%       input files. If not supplied, it tries to guess; the guess could be
%       wrong if you have other files in the path that begin with the same
%       sequence of characters.
%
% To convert to a movie, you might use something like the following (on the
% UNIX command prompt):
%    ffmpeg -r 2 -i framelink%02d.png mymovie.avi
% You will likely need to control the quality, for example with
% "-qscale 5". See the ffmpeg documentation.

% Copyright 2011 by Timothy E. Holy

  if (nargin < 4)
    options = struct;
  end
  options = default(options,...
    'pad_with_zero',false,...
    'overwrite_automatically',false);
  if ~isfield(options,'extension')
    d = dir([basenamein '*']);
    if isempty(d)
      error('Cannot find any input files');
    end
    [~,~,options.extension] = fileparts(d(1).name);
  end
  if options.extension(1) ~= '.'
    options.extension = ['.' options.extension];
  end
  % Check to see if any output files already exist
  d = dir([basenameout '*']);
  if ~isempty(d)
    if options.overwrite_automatically
      delete([basenameout '*']);
    else
      resp = input('Some files with that basename already exist. Delete them? (y/n) ','s');
      if lower(resp(1)) ~= 'y'
        return
      end
    end
  end
  
  n_digits_out = ceil(log10(length(sequence)));
  if options.pad_with_zero
    inputpad = ['%0' num2str(options.pad_with_zero) 'd'];
  else
    inputpad = '%d';
  end
  fmtstr = ['ln -s ' basenamein inputpad options.extension ' ' basenameout '%0' num2str(n_digits_out) 'd' options.extension];
  for i = 0:length(sequence)-1
    [status,msg] = system(sprintf(fmtstr,sequence(i+1),i));
    if (status ~= 0)
      error(msg)
    end
  end
    
