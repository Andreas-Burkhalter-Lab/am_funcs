function mexpp(varargin)
% mexpp: a mex wrapper for nicer development in C++
%
% mexpp tries to ease some of the compiling and linking issues when
% developing more sophisticated MEX files in C++. In conjunction with
% appropriate #include directives within your C++ code, the primary
% features are support for multithreading, support for the high-performance
% C++ matrix library Eigen, and improvements in the MEX API (particularly
% for writing generic templated MEX files and import/export of .mat files).
%
% Syntax/examples:
%
%   mexpp -boostthread myfunction.cpp
%   mexpp -thread myfunction.cpp
% Provides support for a cross-platform threading library, Boost, so that
% you can more easily write portable multithreaded MEX files. The second
% syntax is a shortcut for the first (at least until C++11 threading
% support becomes the default...)
%
%   mexpp -boost myfunction.cpp
% Allow usage of header-only components of Boost, but don't use Boost
% threads.
%
%   mexpp -eigen myfunction.cpp
% Add the #include path for Eigen during compilation
%
%   mexpp -boostthread mythreadpoolserverfunction.cpp
%   mexpp -boostthread -ctrlc mythreadpoolserverfunction.cpp
% If mythreadpoolserverfunction.cpp uses BoostThreadPoolServer, it allows 
% display of progress and the ability to interrupt the MEX file in the
% middle of execution.  The first version uses a dialog box, the second
% allows interrupts via CTRL-C.
%
%   mexpp -standalone ...
% Instead of creating a MEX file, build a standalone application. It can be
% easier to debug (e.g., using a debugger), check for memory leaks (e.g.,
% using valgrind), and profile (again with valgrind) a standalone
% application than a MEX file. You must include a "main" function among
% your source files, which can either be a separate file or can be wrapped
% in "#ifdef MAIN ... #endif" (-standalone defines the symbol MAIN so you
% can test for it).
% -standalone turns on the debugging symbols (-g), which by default turns
% off optimization. Add "-O" to the list of arguments if you want to
% optimize.
%
%   mexpp -setup
% Call this to initialize various search paths. This will happen
% automatically the first time you use mexpp, but if you add new libraries
% you can call this to reconfigure.

% Copyright 2011-2012 by Timothy E. Holy

% Note: Boost.Thread support is provided by "statically-linking" to the
% libary (really, just compiling against the .cpp files). The rationale is
% as follows:
% 1. Matlab itself uses Boost.Thread as a shared library, and chances
%    are it's probably a different version from the one you have on your
%    harddisk. You'll get crashes if you use incompatible versions.
% 2. Even if you carefully match the version you compile with to the
%    version used by your copy of Matlab, someday you might upgrade your
%    copy of Matlab, or you might distribute your MEX file to someone
%    running a different version of Matlab. Since different versions of
%    Matlab use different versions of Boost.Thread, there is no way to
%    ensure that your MEX file will be portable unless you use static
%    linking.
% 3. In my experience, statically linking Boost.Thread increases the
%    size of a MEX file by only a percent or so, which (given the
%    benefits) is negligible.


  %% Parse arguments and initialize
  [options,files] = parse_args(varargin{:});
  % Determine the location of the settings file
  settingsfile = 'mexpp_setup.mat';
  [userpathstr,found] = find_settings_file(settingsfile);
  if ~found
    h = msgbox('This appears to be the first time you''ve run mexpp. You will need to tell mexpp where to find various add-on libraries so they can be added to the #include search paths. If you don''t need to support a particular library, you can click Cancel when it asks you for its location. If you change your mind later, you can always run "mexpp -setup" to incorporate new libraries.','Welcome');
    uiwait(h)
    options.setup = true;
    pause(1); % fixes a bug when working over a slow internet connection
  end
      
  %% Handle the -setup option
  if options.setup
    % Get the tools directory
    toolsdefault = fileparts(which('mexpp'));  % it should be the same one with mexpp in it...
    if exist([toolsdefault filesep 'MatlabTraits.h'],'file')
      mexpptoolsdir = toolsdefault;
    else
      mexpptoolsdir = uigetdir(toolsdefault,'Please find the directory with the MEXppTools .h files (MatlabTraits.h, MatlabIO.h, etc.) (required)');
      if isequal(mexpptoolsdir,0)
        error('The MEXppTools are required (and are distributed with this function, so you probably have them); quitting');
      end
    end
    % Check to see that at least one of the header files is in the chosen
    % directory
    if ~exist([mexpptoolsdir filesep 'MatlabTraits.h'],'file')
      error('The chosen directory does not contain MatlabTraits.h, there must be some error. Quitting.');
    end
    % Get the Eigen directory
    eigendir = get_dir('Eigen','Eigen','EIGENDIR');
    % Get the Boost directory
    boostdir = get_dir('Boost','boost','BOOSTDIR');
    % Get the preferred matopts file for building standalone functions
    if isunix
      matoptsfile = [matlabroot '/bin/matopts.sh'];
      [filename,matoptsdir] = uigetfile('*.sh','Choose the matopts file for compiling standalone executables',matoptsfile);
    else
      matoptsdir = [matlabroot '\bin\win32\mexopts\'];
      % On windows we should determine the right compiler to use
      [~,str] = system('mex -v');
      matchstart = regexp(str,'->    CC','once');
      compilerstart = matchstart + regexp(str(matchstart:end),'= ','once') + 1;
      compilerend = compilerstart+regexp(str(compilerstart:end),'->','once') - 3;
      compilerstr = str(compilerstart:compilerend);
      matoptsfile = [matoptsdir compilerstr 'engmatopts.bat'];
      if ~exist(matoptsfile,'file')
        error('Problem determining the default matopts file, please report bug');
      end
      [filename,matoptsdir] = uigetfile('*engmatopts.bat','Choose the matopts file for compiling standalone executables',matoptsfile);
    end      
    if ~isequal(filename,0)
      matoptsfile = [matoptsdir filename];
    end
    % Save the configuration
    savefile = [userpathstr filesep settingsfile];
    fprintf('Saving settings to %s\n',savefile);
    save(savefile,'mexpptoolsdir','eigendir','boostdir','matoptsfile');
    % On UNIX, check the rpath settings in matopts.sh; if necessary, warn
    % the user about linking with standalone applications
    if isunix
      % Read the matopts.sh file and split it into the different
      % architectures
      fid = fopen(matoptsfile,'r');
      txt = fread(fid,[1 Inf],'*char');
      fclose(fid);
      parts = regexp(txt,';;','split');
      parts = parts(2:end-1);  % first and last parts are not an architecture
      archstr = computer('arch');  % determine this machine's architecture
      for i = 1:length(parts)
        if ~isempty(regexp(parts{i},archstr, 'once'))
          % We've found the right section, now search for "rpath-link"
          if ~isempty(regexp(parts{i},'rpath-link', 'once'))
            warndlg(sprintf('Your %s file uses rpath-link. For standalone executables, this will require you to set the runtime library path, e.g., via LD_LIBRARY_PATH. You can avoid this by editing %s and changing the RPATH variable in the appropriate architecture section (%s) so that it uses -rpath rather than -rpath-link.',matoptsfile,matoptsfile,archstr));
            pause(1)
          end
          break
        end
      end
    end
  end
  
  if isempty(files)
    return
  end
  
  load(settingsfile)
  
  %% Put the tools directory on the search path
  options.flags{end+1} = ['-I' mexpptoolsdir];
  
  %% Convert the options settings into mex flags
  if options.eigen
    options.flags{end+1} = ['-I' eigendir];
  end
  if options.boost
    options.flags{end+1} = ['-I' boostdir];
  end
  if options.boostthread
    % Add the .cpp files we need to compile for thread support
    srcdir = [boostdir filesep 'libs' filesep 'thread' filesep 'src' filesep];
    if isunix
      srcdir = [srcdir 'pthread'];
    elseif ispc
      srcdir = [srcdir 'win32'];
    else
      error('Platform not supported');
    end
    srcdir = [srcdir filesep];
    files = [files {[srcdir 'thread.cpp'],[srcdir 'once.cpp']}];
  end
  if options.ctrlc
    options.flags = [options.flags {'-lut','-DCTRLC'}];
  end
  if options.standalone
    options.flags = [{'-f', matoptsfile,'-DMAIN'} options.flags];
  end
  
  %% Run the compilation
  fprintf('Working...');
  try
    mex(options.flags{:},files{:});
  catch MEXEXCEPT
    fprintf('An error was detected. Possible causes:\n1. A problem in one of your source files\n2. A missing source file\n3. You''re using a library but did not specify the needed flags, e.g., -boost or -lut.\n');
    rethrow(MEXEXCEPT)
  end
  fprintf('Success!\n');
end

function [options,files] = parse_args(varargin)
  options = struct('eigen',false,...
    'boost',false,...
    'boostthread',false,...
    'ctrlc',false,...
    'standalone',false,...
    'setup',false,...
    'flags',{{}});
  files = {};
  for i = 1:length(varargin)
    if ~isempty(varargin{i})
      if varargin{i}(1) == '-'
        % This is an option, so let's try to parse it
        switch lower(varargin{i})
          case '-eigen'
            options.eigen = true;
          case '-boost'
            options.boost = true;
          case {'-boostthread','-thread'}
            options.boost = true;
            options.boostthread = true;
          case '-ctrlc'
            options.ctrlc = true;
          case '-standalone'
            options.standalone = true;
          case '-setup'
            options.setup = true;
          otherwise
            % Can't tell what the option is, but we'd better keep it and
            % pass it on to mex
            options.flags{end+1} = varargin{i};
        end
      else
        % It doesn't begin with '-', so it's presumably a filename
        files{end+1} = varargin{i}; %#ok<AGROW>
      end
    end
  end
end

function [userpathstr,found] = find_settings_file(filename)
  found = exist(filename,'file');
  if found
    userpathstr = which(filename);
    userpathstr = fileparts(userpathstr);
  else
    userpathstr = userpath;
    % On UNIX systems, there can be multiple entries, e.g., if the user has a
    % directory named "matlab" in the home folder
    if isunix
      userpathc = regexp(userpathstr,':','split');
      userpathstr = userpathc{1};
    end
  end
end

function dirname = get_dir(showname,checkname,envname)
  % Get the directory
  if (nargin > 2)
    dirname = getenv(envname);
    if ~isempty(dirname) && exist([dirname filesep checkname],'dir')
      answer = questdlg(sprintf('I found the environment variable %s with a suitable setting (%s). Use this for the %s search path?',envname,dirname,showname),'Environment variable found','Yes','No','Yes');
      pause(1)
      if strcmpi(answer,'yes')
        return
      end
    end
  end
  while true
    dirname = uigetdir('',sprintf('If you have a download of %s, please help me find it (hint: it should be a directory containing a subdirectory named "%s")',showname,checkname));
    if isequal(dirname,0)
      dirname = '';
      break
    else
      % Check for the required subdirectory
      if ~exist([dirname filesep checkname],'dir')
        h = msgbox(sprintf('There is no sub-folder named "%s" in this directory. Please try again, or click Cancel if you don''t use %s.',checkname,showname),'Not found','warn');
        uiwait(h)
        pause(1)
      else
        % It was found successfully
        break
      end
    end
  end
end
