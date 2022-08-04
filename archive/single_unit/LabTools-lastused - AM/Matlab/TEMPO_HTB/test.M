function htbtest(filename, db)
%HTBTEST - Sample program to exercise htb...() functions
%
% This program uses some of htb functions to display a database
% in an HTB file.  It opens the TEST.HTB file and displays the
% information about each database.
%
% It copies a database to TEST1.HTB.
%
% It also displays the data of an epoch with htbShow().
%
% The following functions are demonstrated:
%   htbOpen(), htbCount(), htbGetHd(), htbGetEp(), htbClose()
%   htbOpenW(), htbWrite(), htbGetDa()
%   htbShow() is called to display an epoch from a database.

filename = 'test.htb';              % Let's look at this file

fid = htbOpen(filename);            % Open our HTB file
ndbs = htbCount(fid);               % Get number of databases in it
t = sprintf('File %s has %d databases', filename, ndbs);
disp(t);

for database = 1:ndbs               % Loop though each one...
    fprintf('  Reading database header %d: ', database);
    htb{database} = htbGetHd(fid, database);   % HTB header
    h = htb{database};              % for convenience
    fprintf('%s (%s)\n', h.title, h.pro_file);
    hertz = h.speed_units / h.speed;    % see p366 in TEMPO v9 manual
    binwidth = (h.skip + 1) / hertz;
    epoch_start = (-h.offset) * binwidth;
    epoch_period = h.period * binwidth;
    epoch_extent = h.extension * binwidth;
    fprintf('    Binwidth=%0.3f(sec),  Window=%0.3f(sec)-%0.3f(sec)   Extent=%0.3f(sec)\n', ...
        binwidth, epoch_start, epoch_start + epoch_period, epoch_extent);
    fprintf('    Reading database %d (%s) epoch 1..\n', database, htb{database}.title);
    epoch{database} = htbGetEp(htb{database}, 1); % Get first epoch
    end

e1 = htbGetDa(htb{1});              % Get all epoch data for database 1

err = htbClose(fid);                % Close HTB file

% SHOW HOW TO WRITE A DATABASE TO ANOTHER FILE.
% THE FILE IS CREATED, IF NECESSARY.
% THE HEADER AND DATABASE ARE APPENDED TO THE FILE.
% DATABASES ALREADY IN THE FILE ARE LEFT UNDISTURBED.

fprintf('\nCopying the first database in %s to TEST1.HTB..', filename);

fid1 = htbOpenW('test1.htb');       % Open another file
err = htbWrite(fid1, htb{1}, e1);   % Write header & database to another file
err = htbClose(fid1);               % Close second file

%whos                               % display variables

htb{1}                              % Display a database header for fun

% GET A LITTLE FANCIER AND DISPLAY SOME DATA

htbshow(filename, 1, 1);            % Display db 1 epoch 1
