function fmb_save = fmb_handscoring_log(in_fmbs, cont)
% function fmb_handscoring_log allows user to blindly analyze fmb data
%
% Syntax: first arg 'in_fmbs' should either:
%         A)   be the char array 'new' for an input dialog
%         B)   point to a .fmbs (fmb save) file 
%         C)   be a compatible cell array of filename strings (full path), 
%              with the first in each cell referring to the movie 
%              ('.yuv') file and the second referring to the '.fmbn' file
%         D)   the char array 'key' to return to screen the decoding of a
%              completed trial
%
% Copyright 2008 Julian P Meeks (Timothy Holy Laboratory)
%
% Version History
% 2008_07_25: wrote it (JPM)

%% lets start by checking inputs and initializing in_fmb if 'new'
%
if ~iscell(in_fmbs)
    if ~isempty(strmatch(in_fmbs, 'new'))
        movie_filenames = []; fmb_filenames = [];
        while size(movie_filenames,2) ~= size(fmb_filenames,2) ... 
            || isempty(movie_filenames) ...
            || isempty(fmb_filenames)
        movie_filenames = UIGetFiles('*.yuv', 'Please select movie files (in order)', pwd);
        fmb_filenames = UIGetFiles('*.fmbn', 'Please select .fmbn files (in order)', pwd);
        end
        for idx = 1:size(movie_filenames,2)
            in_fmb{idx}(1) = movie_filenames(idx);
            in_fmb{idx}(2) = fmb_filenames(idx);
        end
        cont = false;
    elseif exist(in_fmbs, 'file')
        load(in_fmbs, '-mat');
        if ~exist('fmb_save', 'var')
            error('filename does not contain a compatible fmb_save variable');
        else
            in_fmb = fmb_save.fmb_files;
            if isfield(fmb_save, 'completed')
                cont = true;
                fmb_complete = fmb_save.completed;
            end
        end
    end
else
    cont = FALSE;
    in_size = size(in_fmb, 2);
    for idx = 1:in_size
        if ~iscellstr(in_fmbs{idx})
            error('first argument must be a cell array of cell array of strings.');
        elseif size(in_fmbs{idx})~=2
            error('cell array of strings does not contain paired filenames.');
        end
    end
    in_fmb = in_fmbs;
end

%% Now prepare filenames
infmb_size = size(in_fmb,2);
for idx = 1:infmb_size
    movie_filenames{idx} = in_fmb{idx}{1};
    fmb_filenames{idx} = in_fmb{idx}{2};
end

%% if new, set up save struct
if ~exist('fmb_save', 'var')
    fmb_save = struct;
    fmb_save.fmb_files = in_fmb;
end

%% start up the loop
if ~isfield(fmb_save, 'completed')
    fmb_save.completed = [];
end

passload = false;

while size(fmb_save.completed,2) < infmb_size
    
    if exist('answer', 'var')
        if isempty(strmatch(answer, 'No'))
            % Choose a file to use @ random
            nextfile = ceil(abs(rand)*infmb_size);
            if cont
                while ~isempty(intersect(nextfile, fmb_save.completed))
                    nextfile = ceil(abs(rand)*infmb_size);
                end
            end
        else
            passload = true;
        end
    else
        nextfile = ceil(abs(rand)*infmb_size);
        if cont
            while ~isempty(intersect(nextfile, fmb_complete))
                nextfile = ceil(abs(rand)*infmb_size);
            end
        end
    end
    % ASSERT: the 'nextfile' index is a not-previously-completed file

    
    % load movie file
    % note! this needs to be fixed to automatically determine the movie frame
    % size!! otherwise it will barf if a strange file size is used!
   if passload == false
    % load fmbn file
     load(in_fmb{nextfile}{2}, '-mat');
     if ~exist('fmb', 'var')
         error([in_fmb{nextfile}{2} ' does not contain a proper ''fmb'' variable.']);
     end
    % load movie file
       if exist('mov', 'var')
           clear mov;
       end
       %mov = yuv2mov(in_fmb{nextfile}{1}, 240, 352, '420');
   end
    % run fmb_moviewatch
    try
        fmb_moviewatch(mov, fmb);
    catch
        lerr = lasterror;
        if strmatch(lerr.identifier, 'MATLAB:UndefinedFunction')
            delete(figure(1));
            delete(figure(2));
        else
            return;
        end
    end

    % Ask if completed:
    answer = questdlg('Did you successfuly complete the analysis?',...
                  'Analysis complete?',...
                  'Yes', 'No', 'Yes');
    if ~isempty(strmatch(answer, 'Yes'))
        previous = [];
        passload = false;
        fmb_save.handscore_categories = {'Sniffs Only', 'Sniff->Mounts', 'Licks/bites only', 'Licks/bites-mounts', 'Touch/mounts only'};
        input_values = []
        for idx = 1:4 
            while isempty(input_values)
                if isempty(previous)
                    input_values{idx} = inputdlg(fmb_save.handscore_categories, ['Controller #' num2str(idx)]);
                else
                    input_values{idx} = inputdlg(fmb_save.handscore_categories, ['Controller #' num2str(idx)], 1, previous);
                end
                for in_idx = 1:size(input_values{idx},1)
                    if isempty(str2num(input_values{idx}{in_idx}))
                        previous = input_values{idx};
                        input_values = [];
                        break;
                    else
                        integer_values(idx, in_idx) = str2num(input_values{idx}{in_idx});
                    end
                end
            end
            previous = [];
            input_values = [];
        end

        if cont
            fmb_save.handscores{end+1} = integer_values;
        else
            fmb_save.handscores{1} = integer_values;
        end

        if cont
            fmb_save.completed(end+1) = nextfile;
        else
            fmb_save.completed(1) = nextfile;
            cont = true;
        end
        save([datestr(now, 'yyyymmdd') 'inprogress' '.fmbs'], 'fmb_save', '-mat');
    end
end

save([datestr(now, 'yyyymmdd') datestr(now,'HHss') '.fmbs'], 'fmb_save', '-mat');
delete([datestr(now, 'yyyymmdd') 'inprogress' '.fmbs']);

end