% Celf6 Juveniles 2013 - Krystal's Recordings Folder

% So, the first thing, is I'm going to make sonograms with no threshold
% in order to properly pick threshold using freqrange [50000 100000]

make_sonograms; % plot.threshold=0; plot.freqrange =[50000 100000]

% Code for setting manual thresholds:

% The newer problem: Krystal didn't have all these files coded
% consistently. She says that files 003 - 006 are all HAB for the same
% animal. Let's see what all the lengths are first:
wavnames = dir('*.WAV');
sizes = zeros(1,numel(wavnames));
for i = 1:numel(wavnames)
    s = wavread(wavnames(i).name,'size');
    sizes(i) = s(1);
    clear s;
end
sizes = sizes./250000;


% The problem: With all the background noise from digging in bedding and
% running around in the cage, even if low frequency (<10kHz), setting a
% high threshold will probably filter out my whistles as anything in the
% SNG which is less than threshold will be set to zero

% So instead, sample sonograms with no threshold, on 50-10 kHz, then set
% threshold based on the time vs. max amplitude graph as follows:


% Manual thresholding
% make all the sonograms first using:
% sngparms.plot = 0;
% sngparms.threshold = 0;
% sngparms.nfreq = 256;
% sngparms.freqrange = [50000 100000];
% Then process them by running the script below:

wavnames = dir('*.WAV');
sngnames = dir('*.SNG');

bins = cell(1,numel(sngnames));

for i = 1:numel(wavnames)
s = wavread(wavnames(i).name,'size');
s = s(1);
s = s/250000;
binn = s/1;
bins{i} = linspace(1,s,binn);
bins{i} = [bins{i}(1:end-1);bins{i}(2:end)];
end
x = cell(1,10);
y = cell(1,numel(sngnames));

for f = 1:numel(sngnames);
% by time
% This finds the max amplitude x time in the sonogram
% in each time bin over 20 kHz to 120 kHz
if length(bins{f}) > 40
for i = 20:39 % 20 second time window
    [sng,header] = ReadSonogram(sngnames(f).name,[bins{f}(1,i) bins{f}(2,i)]);
    sng = abs(sng);
    m = max(sng);
    x{i-19} = m(find(full(m)));
    clear m;
    clear sng;
end

y{f} = cat(2,x{:});
clear x;
else
    y{f} = 0;
end
end

% forplot has the appearance of mxn where m = animal, n is max
% this is the max amp in all rows in every msec column over the 20 sec interval

% Now time to manually enter the threshholds
% Initialize the plots structure

plots(numel(sngnames)).id = wavnames(end).name;
plots(numel(sngnames)).threshold = 0;

% Manual thresholding interface:
for i = 1:numel(sngnames)
    test = 0;
    if ~isempty(find(y{i})) == 1
    while test == 0
        close;
        h = figure;
        % time axis
        k = find(y{i});
        p = semilogy(k,y{i});
        % Make the graph interactive for grabbing data:
        dcm_obj = datacursormode(h);
        set(dcm_obj,'DisplayStyle','window','Enable','on');
        % plot max amp 70 - 100 kHz, on bins{i}(1,20) - bins{i}(2,29)
        
        plots(i).id = wavnames(i).name;
        
        % On plot, use data cursor to estimate the threshold
        disp('Okiedoke, click on where you think the threshold is. Then press Enter.');
        pause; % Now wait for the user. You have to make sure to go back to the command window to hit Enter!
        % If you hit Enter in the figure window, it saves the threshold you
        % clicked and draws the next graph. There should be a way to fix
        % this.
        c_info = getCursorInfo(dcm_obj);
        
        plots(i).threshold = c_info.Position(2); 
        
        % Now, draw a red line and ask the user if it's ok.
        hold on;
        l = line([k(1) k(end)],[plots(i).threshold plots(i).threshold],'Color','red');
        
        ask = input('Does it look good? Enter 0 for no and 1 for yes: ');
        
        while ask == 0
            w = input('Well, do you want: (1) Move the line based on input, (2) Redraw and click around again?  ');
            if w == 1
                plots(i).threshold = input('Enter desired threshold: ');
                l = line([k(1) k(end)],[plots(i).threshold plots(i).threshold],'Color','red');
                ask = input('(1) Look good? (0) or not? ');
            elseif w == 2
                disp('Okiedoke, click on where you think the threshold is. Then press Enter.');
                pause; % Now wait for the user.
                c_info = getCursorInfo(dcm_obj);
                plots(i).threshold = c_info.Position(2);
                l = line([k(1) k(end)],[plots(i).threshold plots(i).threshold],'Color','red');
                ask = input('(1) Look good? (0) or not? ');
            end
        end
        
        if ask == 1
            test = 1;
        end
        
        
    end
    close;
    else
        plots(i).id = wavnames(i).name;
        plots(i).threshold = {'no data'};
    end
    
end

save('manualthresholds.mat','plots');

% 8/11/2013
% So I took a look at all the sonograms. The whistles are clear, but a
% little patchy - meaning we may still have a threshold that's a little to
% high.
% Tweaking file 020.WAV. First by 0.5*plots(i).threshold
% If that works, then go back and remake all the sonograms with half the
% threshold.

% Try whistimes with 10 min sonogram clip, otherwise it'll crashb:
[sng,header]=ReadSonogram('C6JUV_020_500.WAV.sng');
f = linspace(0,header.scanrate/2,header.nfreq);
t = linspace(0,header.nscans/header.scanrate,header.columnTotal);
dt = t(end)/header.columnTotal; % time step in seconds per column
% Create a 10 minute sonogram
ncols = round(600/dt); % 600 seconds / dt seconds per column = x columns
ds = header.scanrate; % number of scans per second
nscans = header.scanrate*600; % number of scans in 10 minutes
newheader = header;
newheader.nscans = nscans;
newheader.columnTotal = ncols;
newheader.tacq = 600;
sng = sng(:,1:ncols);
twhis = whistimes(sng,newheader,whistimesdefaults);

% Ok, processing in 10 minute blocks seems to work. Still, the thresholding
% sucks. Let's make a different one. try 600?
sngparms = struct('plot',0,'threshold',600,...
                   'freqrange',[40000 120000],'nfreq',256);



% after trying a bunch, I think I have to go with 600. 

% Now, let's calculate twhis for each 10 minute bin.

sngnames = dir('*.SNG');
sngs(numel(sngnames)).sngname = sngnames(numel(sngnames)).name;
sngs(numel(sngnames)).id = [];
sngs(numel(sngnames)).filetype = {};
sngs(numel(sngnames)).tbins = [];
sngs(numel(sngnames)).twhis1 = [];
sngs(numel(sngnames)).twhis2 = [];
sngs(numel(sngnames)).twhis3 = [];
sngs(numel(sngnames)).dt1 = [];
sngs(numel(sngnames)).dt2 = [];
sngs(numel(sngnames)).dt3 = [];
sngs(numel(sngnames)).totvtime1 = [];
sngs(numel(sngnames)).totvtime2 = [];
sngs(numel(sngnames)).totvtime3 = [];
sngs(numel(sngnames)).scores1 = {};
sngs(numel(sngnames)).scores2 = {};
sngs(numel(sngnames)).scores3 = {};
sngs(numel(sngnames)).keys1 = {};
sngs(numel(sngnames)).keyscores1 = {};
sngs(numel(sngnames)).keydt1 = {};
sngs(numel(sngnames)).totkeytime1 = {};
sngs(numel(sngnames)).keys2 = {};
sngs(numel(sngnames)).keyscores2 = {};
sngs(numel(sngnames)).keydt2 = {};
sngs(numel(sngnames)).totkeytime2 = {};
sngs(numel(sngnames)).keys3 = {};
sngs(numel(sngnames)).keyscores3 = {};
sngs(numel(sngnames)).keydt3 = {};
sngs(numel(sngnames)).totkeytime3 = {};



% TEST FILE WHISTIMES
for i = [2,7,9,11,13,15,17,19]
[sng,header]=ReadSonogram(sngnames(i).name);
f = linspace(0,header.scanrate/2,header.nfreq);
t = linspace(0,header.nscans/header.scanrate,header.columnTotal);
dt = t(end)/header.columnTotal;
ds = header.scanrate;
ncols = round(600/dt);
% Now, grab indices for 3 10 minute bins.
bins = [1,(ncols + 1), (2*ncols + 1);...
        ncols, (2*ncols), size(sng,2)];
totcols = (bins(2,:) - bins(1,:))+1;
% Now, find the time correspondences for each 10 minute bin so we
% can convert the twhis later
tbins = [dt, (ncols+1)*dt, (2*ncols+1)*dt;...
         ncols*dt, (2*ncols)*dt, header.columnTotal*dt];
tott = tbins(2,:) - tbins(1,:);

header1 = header;
header1.nscans = round(ds*tott(1));
header1.columnTotal = totcols(1);
header1.tacq = tott(1);
sng1 = sng(:,bins(1,1):bins(2,1));

header2 = header;
header2.nscans = round(ds*tott(2));
header2.columnTotal = totcols(2);
header2.tacq = tott(2);
sng2 = sng(:,bins(1,2):bins(2,2));

header3 = header;
header3.nscans = round(ds*tott(3));
header3.columnTotal = totcols(3);
header3.tacq = tott(3);
sng3 = sng(:,bins(1,3):bins(2,3));

twhis1 = whistimes(sng1,header1,whistimesdefaults);
twhis2 = whistimes(sng2,header2,whistimesdefaults);
twhis3 = whistimes(sng3,header3,whistimesdefaults);

sngs(i).sngname = sngnames(i).name;
sngs(i).twhis1 = twhis1;
sngs(i).twhis2 = twhis2;
sngs(i).twhis3 = twhis3;
sngs(i).tbins = tbins;
sngs(i).filetype = 'test';


clearvars -except sngnames sngs

end
% Ok, let's 'believe' it for now, even though I can see there are some
% mistakes. Let's try to sync the keys.

% Calculate the twhis for the hab files, for completeness. Though they
% should all be 'zero' basically, and there's only 1 bin.
% HAB FILE WHISTIMES
for i = [1,3,4,5,6,8,10,12,14,16,18]
[sng,header]=ReadSonogram(sngnames(i).name);
twhis1 = whistimes(sng,header,whistimesdefaults);
sngs(i).sngname = sngnames(i).name;
sngs(i).twhis1 = twhis1;
sngs(i).filetype = 'hab';
end


% Now let's import, clean up the logfiles and try to do some syncing of
% audio to the scored behaviors.

% We'll make the keystroke key later. First 
% 1. Read the CSVs
% 2. Trim them? Load logfile into the sngs structure
% 3. Grab for each keystroke the start and end times. Make a
% sngs(i).keystroke field for each keystroke

% Before I do it, ran make_sonograms again! UGH - the weird lines
% throughout the sonogram are causing problems. Set freqrange = [50000 120000]
% with custom thresholds again.

% 8/13/2013 This worked out really well! The sonograms look a lot cleaner
% than before. And the whistles pop out nicely.

% Let's practice loading and editing CSVs

cd Celf6' JuvPlay log files'\;
a = dir('*.CSV');
csvfile_id = zeros(1,numel(a));
csvfile_bin = zeros(1,numel(a));

% Now let's make a vector of the ids and bin of every file
%
for i = 1:numel(a);
    n=a(i).name;
    n=strrep(n,' log 1.csv','');
    n=strrep(n,'Scored_JuvPlay_C6-','');
    % now it's of the form ID_BIN;
    n=cat(2,n);
    k=find(n=='_');
    csvfile_id(i) = str2num(n(1:(k-1)));
    csvfile_bin(i) = str2num(n((k+1):end));
end
  
clear n i k;
% Great, now we have three variables :
% a :  the filenames
% csvfile_id : the file ids - match to sngs(i).id
% csvfile_bin: the 10 minute time bin - matches to twhis1, twhis2, twhis3


% Now manually edit the sngs structure to add ids from Krystal's
% spreadsheet

% I imported all the ids for the 19 files into x
x = [1,1,6,6,6,6,6,10,10,9,9,6,6,14,14,13,13,12,12];
for i = 1:numel(sngs)
    sngs(i).id = x(i);
end
clearvars -except a csv* sng*

% Ok, now we're ready to start getting to work.

% Real code for reading the csvs starts here:

% READ THE CSVS
for i = 1:numel(a)
x{i} = textread(a(i).name,'%s','delimiter','\n');
% This makes a cell array, where each entry is each line from the csvfile
end
csvs = cell(1,length(x));
for j = 1:length(x)
y = x{j}';
row = cell(1,length(y));
    for i = 1:length(y)
        row{i} = regexp(y{i},',','split');     
    end
d = cat(1,row{3:end});
csvs{j} = d;
clear y d;
end

clearvars -except a csv* sng*

% Now, let's add the csvs to the sngs structure. First, let's find all the
% ids:

x = cat(2,sngs(:).id);


% ADD THE CSVS TO THE SNGS STRUCTURE
for i = 1:length(x)

    if strcmp(sngs(i).filetype,'test') == 1 
        % only write if it's a test file
    
        binidx1 = find(csvfile_id == x(i) & csvfile_bin == 1);
        binidx2 = find(csvfile_id == x(i) & csvfile_bin == 2);
        binidx3 = find(csvfile_id == x(i) & csvfile_bin == 3);
        
        sngs(i).scores1 = csvs{binidx1};
        sngs(i).scores2 = csvs{binidx2};
        sngs(i).scores3 = csvs{binidx3};
        
        clear binidx*
        
    
    end
    


end

% Great that totally worked!!

clearvars -except sngs;

% So now sngs is loaded with everything

% Now the analysis and sorting.

% First Time Bin
for i = 1:numel(sngs)
     if strcmp(sngs(i).filetype,'test') == 1 
        keys = sngs(i).scores1(:,1);
        startstop = sngs(i).scores1(:,2);
        times = sngs(i).scores1(:,3);

        keys = keys';

        startstop = startstop';

        times = times';
            for j = 1:length(times);
              timesdouble(j) = str2num(times{j});
            end
        clear times;



        % Trim 'begins' and 'ends'

        keys = keys(2:(end-1));
        startstop = startstop(2:(end-1));
        timesdouble = timesdouble(2:(end-1));

        [keycodes,keynames] = grp2idx(keys);
        [startid, startname] = grp2idx(startstop);
        keyscores = cell(1,length(keynames));

        for j = 1:length(keynames)
    
            k = find(keycodes == j);
            x = startid(k);
            t = timesdouble(k);
            y = [t(find( x == 1));t(find(x == 2))];
            y(3,1:size(y,2)) = j;
            keyscores{j} = y;
            clear k x t y
        end  
    
    
% Great.

% Now add this info to sngs

         sngs(i).keys1 = keynames;
         sngs(i).keyscores1 = keyscores;
         
         clear key* start* times*;
     end
 
end
% Second Time Bin
for i = 1:numel(sngs)
     if strcmp(sngs(i).filetype,'test') == 1 
        keys = sngs(i).scores2(:,1);
        startstop = sngs(i).scores2(:,2);
        times = sngs(i).scores2(:,3);

        keys = keys';

        startstop = startstop';

        times = times';
            for j = 1:length(times);
              timesdouble(j) = str2num(times{j});
            end
        clear times;



        % Trim 'begins' and 'ends'

        keys = keys(2:(end-1));
        startstop = startstop(2:(end-1));
        timesdouble = timesdouble(2:(end-1));

        [keycodes,keynames] = grp2idx(keys);
        [startid, startname] = grp2idx(startstop);
        keyscores = cell(1,length(keynames));

        for j = 1:length(keynames)
    
            k = find(keycodes == j);
            x = startid(k);
            t = timesdouble(k);
            y = [t(find( x == 1));t(find(x == 2))];
            y(3,1:size(y,2)) = j;
            keyscores{j} = y;
            clear k x t y
        end  
    
    
% Great.

% Now add this info to sngs

         sngs(i).keys2 = keynames;
         sngs(i).keyscores2 = keyscores;
         
         clear key* start* times*;
     end
 
end
% Third Time Bin
for i = 1:numel(sngs)
     if strcmp(sngs(i).filetype,'test') == 1 
        keys = sngs(i).scores3(:,1);
        startstop = sngs(i).scores3(:,2);
        times = sngs(i).scores3(:,3);

        keys = keys';

        startstop = startstop';

        times = times';
            for j = 1:length(times);
              timesdouble(j) = str2num(times{j});
            end
        clear times;



        % Trim 'begins' and 'ends'

        keys = keys(2:(end-1));
        startstop = startstop(2:(end-1));
        timesdouble = timesdouble(2:(end-1));

        [keycodes,keynames] = grp2idx(keys);
        [startid, startname] = grp2idx(startstop);
        keyscores = cell(1,length(keynames));

        for j = 1:length(keynames)
    
            k = find(keycodes == j);
            x = startid(k);
            t = timesdouble(k);
            y = [t(find( x == 1));t(find(x == 2))];
            y(3,1:size(y,2)) = j;
            keyscores{j} = y;
            clear k x t y
        end  
    
    
% Great.

% Now add this info to sngs

         sngs(i).keys3 = keynames;
         sngs(i).keyscores3 = keyscores;
         
         clear key* start* times*;
     end
 
end

clearvars -except sngs;

% Scoring voc vs. key sync

% 1:        twhis start or end within keyscore start&end OVERLAP
% 2:        twhis end within 1 sec of keyscore start     PRECEDE
% 3:        twhis start within 1 sec of keyscore end     FOLLOW
% 0:        none of the above

% SCORE THE VOC BEHAVIOR SYNC <--- OLD, REWRITE



% Not the best scoring scheme. First some calculations

% Calculate whistle duration times

for i = 1:numel(sngs)
    
  if ~isempty(sngs(i).twhis1)
    sngs(i).dt1 = sngs(i).twhis1(2,:) - sngs(i).twhis1(1,:);
    sngs(i).totvtime1 = sum(sngs(i).dt1);
  end
  if ~isempty(sngs(i).twhis2)
    sngs(i).dt2 = sngs(i).twhis2(2,:) - sngs(i).twhis2(1,:);
    sngs(i).totvtime2 = sum(sngs(i).dt2);
  end
  if ~isempty(sngs(i).twhis3)
    sngs(i).dt3 = sngs(i).twhis3(2,:) - sngs(i).twhis3(1,:);
    sngs(i).totvtime3 = sum(sngs(i).dt3);
  end
end


% Calculate behavior duration times

for i = 1:numel(sngs)
   if ~isempty(sngs(i).keys1)
       for j = 1:length(sngs(i).keyscores1)
        sngs(i).keyscores1{j} = sngs(i).keyscores1{j}(1:2,:);
        sngs(i).keydt1{j} = sngs(i).keyscores1{j}(2,:) - sngs(i).keyscores1{j}(1,:);
        sngs(i).totkeytime1{j} = sum(sngs(i).keydt1{j}(:));
       end
   end
   if ~isempty(sngs(i).keys2)
       for j = 1:length(sngs(i).keyscores2)
        sngs(i).keyscores2{j} = sngs(i).keyscores2{j}(1:2,:);
        sngs(i).keydt2{j} = sngs(i).keyscores2{j}(2,:) - sngs(i).keyscores2{j}(1,:);
        sngs(i).totkeytime2{j} = sum(sngs(i).keydt2{j}(:));
       end
   end
   if ~isempty(sngs(i).keys3)
       for j = 1:length(sngs(i).keyscores3)
        sngs(i).keyscores3{j} = sngs(i).keyscores3{j}(1:2,:);
        sngs(i).keydt3{j} = sngs(i).keyscores3{j}(2,:) - sngs(i).keyscores3{j}(1,:);
        sngs(i).totkeytime3{j} = sum(sngs(i).keydt3{j}(:));
       end
   end
end


% SCORING FOR TIME BIN 1

% USING THE NEWLY MADE TESTS STRUCTURE THAT IS 1-6 of just the tests

% FIRST TIME BIN STUFF
for i = 1:numel(tests)
    if ~isempty(tests(i).keys1)
    
        x = tests(i).keyscores1;
        twhis = tests(i).twhis1;
        dt = tests(i).dt1;
        
        for key = 1:length(x);
            y = x{key};
            
            for bout = 1:size(y,2) % loop over each bout for keyth key
            
                bstart = y(1,bout);
                bend = y(2,bout);
                
                
                overlaps = find((twhis(1,:)>= bstart & twhis(2,:) <= bend) | ...
                                (twhis(1,:)<= bstart & twhis(2,:) >= bstart)| ...
                                (twhis(1,:)<= bend & twhis(2,:) >= bend));
                            
                precedes = find((twhis(1,:)<= bstart & twhis(2,:) <= bstart) & ...
                                (abs(twhis(2,:) - bstart) <= 1));
                            
                follows = find((twhis(1,:) >= bend) & (twhis(2,:) >= bend) & ...
                                (abs(twhis(1,:) - bstart) <= 1));
                            
                            
                oversum{key}(bout) = sum(dt(overlaps));
                precedesum{key}(bout) = sum(dt(precedes));
                followsum{key}(bout) = sum(dt(follows));
            
                clear bstart bend overlaps precedes follows
            end
            clear y;
        end
        
        clear x twhis dt;
        
        overs1{i} = oversum;
        precs1{i} = precedesum;
        folls1{i} = followsum;
        
        clear oversum precedesum followsum;
    end
end

% Second time bin
for i = 1:numel(tests)
    if ~isempty(tests(i).keys2)
    
        x = tests(i).keyscores2;
        twhis = tests(i).twhis2;
        dt = tests(i).dt2;
        
        for key = 1:length(x);
            y = x{key};
            
            for bout = 1:size(y,2) % loop over each bout for keyth key
            
                bstart = y(1,bout);
                bend = y(2,bout);
                
                
                overlaps = find((twhis(1,:)>= bstart & twhis(2,:) <= bend) | ...
                                (twhis(1,:)<= bstart & twhis(2,:) >= bstart)| ...
                                (twhis(1,:)<= bend & twhis(2,:) >= bend));
                            
                precedes = find((twhis(1,:)<= bstart & twhis(2,:) <= bstart) & ...
                                (abs(twhis(2,:) - bstart) <= 1));
                            
                follows = find((twhis(1,:) >= bend) & (twhis(2,:) >= bend) & ...
                                (abs(twhis(1,:) - bstart) <= 1));
                            
                            
                oversum{key}(bout) = sum(dt(overlaps));
                precedesum{key}(bout) = sum(dt(precedes));
                followsum{key}(bout) = sum(dt(follows));
            
                clear bstart bend overlaps precedes follows
            end
            clear y;
        end
        
        clear x twhis dt;
        
        overs2{i} = oversum;
        precs2{i} = precedesum;
        folls2{i} = followsum;
        
        clear oversum precedesum followsum;
    end
end
 
% Third time bin
for i = 1:numel(tests)
    if ~isempty(tests(i).keys3)
    
        x = tests(i).keyscores3;
        twhis = tests(i).twhis3;
        dt = tests(i).dt3;
        
        for key = 1:length(x);
            y = x{key};
            
            for bout = 1:size(y,2) % loop over each bout for keyth key
            
                bstart = y(1,bout);
                bend = y(2,bout);
                
                
                overlaps = find((twhis(1,:)>= bstart & twhis(2,:) <= bend) | ...
                                (twhis(1,:)<= bstart & twhis(2,:) >= bstart)| ...
                                (twhis(1,:)<= bend & twhis(2,:) >= bend));
                            
                precedes = find((twhis(1,:)<= bstart & twhis(2,:) <= bstart) & ...
                                (abs(twhis(2,:) - bstart) <= 1));
                            
                follows = find((twhis(1,:) >= bend) & (twhis(2,:) >= bend) & ...
                                (abs(twhis(1,:) - bstart) <= 1));
                            
                            
                oversum{key}(bout) = sum(dt(overlaps));
                precedesum{key}(bout) = sum(dt(precedes));
                followsum{key}(bout) = sum(dt(follows));
            
                clear bstart bend overlaps precedes follows
            end
            clear y;
        end
        
        clear x twhis dt;
        
        overs3{i} = oversum;
        precs3{i} = precedesum;
        folls3{i} = followsum;
        
        clear oversum precedesum followsum;
    end
end

% YAY!!!!!!

% Ok, now the Rate (msec voc/sec behavior) calculations

% First Time Bin
for i = 1:numel(overs1)
    for key = 1:numel(overs1{i})
        overrate1{i}(key) = (sum(overs1{i}{key})*1000)/tests(i).totkeytime1{key};
        precrate1{i}(key) = (sum(precs1{i}{key})*1000)/tests(i).totkeytime1{key};
        follrate1{i}(key) = (sum(folls1{i}{key})*1000)/tests(i).totkeytime1{key};
    end
end

% Second Time Bin
for i = 1:numel(overs2)
    for key = 1:numel(overs2{i})
        overrate2{i}(key) = (sum(overs2{i}{key})*1000)/tests(i).totkeytime2{key};
        precrate2{i}(key) = (sum(precs2{i}{key})*1000)/tests(i).totkeytime2{key};
        follrate2{i}(key) = (sum(folls2{i}{key})*1000)/tests(i).totkeytime2{key};
    end
end

% Third Time Bin
for i = 1:numel(overs3)
    for key = 1:numel(overs3{i})
        overrate3{i}(key) = (sum(overs3{i}{key})*1000)/tests(i).totkeytime3{key};
        precrate3{i}(key) = (sum(precs3{i}{key})*1000)/tests(i).totkeytime3{key};
        follrate3{i}(key) = (sum(folls3{i}{key})*1000)/tests(i).totkeytime3{key};
    end
end
clear i key bout;
% Find the instersection of keys scored in all files in each bin.

for i = 1:6
    keys1{i} = cat(2,tests(i).keys1{:});
    keys2{i} = cat(2,tests(i).keys2{:});
    keys3{i} = cat(2,tests(i).keys3{:});
end


keys1int = intersect(intersect(intersect(intersect(intersect(keys1{1},keys1{2}),keys1{3}),keys1{4}),keys1{5}),keys1{6});
keys2int = intersect(intersect(intersect(intersect(intersect(keys2{1},keys2{2}),keys2{3}),keys2{4}),keys2{5}),keys2{6});
keys3int = intersect(intersect(intersect(intersect(intersect(keys3{1},keys3{2}),keys3{3}),keys3{4}),keys3{5}),keys3{6});


keys1int = keys1{1};
for i = 2:length(keys1)
    keys1int = intersect(keys1int,keys1{i});
end
keys2int = keys2{1};
for i = 2:length(keys2)
    keys2int = intersect(keys2int,keys2{i});
end
keys3int = keys2{1};
for i = 2:length(keys3)
    keys3int = intersect(keys3int,keys3{i});
end


save overprecfoll_withkeystrokeintersections.mat;
% OK great. Now, look at msec/sec activity over, preceding and folliwng
% each keystroke in keys intersections.

% Time Bin 1
for animal = 1:length(tests)

    for key = 1:length(keys1int)
        idx = find(keys1{animal} == keys1int(key));
        
        overrate1_clip{animal}(key) = overrate1{animal}(idx);
        precrate1_clip{animal}(key) = precrate1{animal}(idx);
        follrate1_clip{animal}(key) = follrate1{animal}(idx);
        
        clear idx;
    end
        

end

% Time Bin 2
for animal = 1:length(tests)

    for key = 1:length(keys2int)
        idx = find(keys2{animal} == keys2int(key));
        
        overrate2_clip{animal}(key) = overrate2{animal}(idx);
        precrate2_clip{animal}(key) = precrate2{animal}(idx);
        follrate2_clip{animal}(key) = follrate2{animal}(idx);
        
        clear idx;
    end
        

end

% Time Bin 3
for animal = 1:length(tests)

    for key = 1:length(keys3int)
        idx = find(keys3{animal} == keys3int(key));
        
        overrate3_clip{animal}(key) = overrate3{animal}(idx);
        precrate3_clip{animal}(key) = precrate3{animal}(idx);
        follrate3_clip{animal}(key) = follrate3{animal}(idx);
        
        clear idx;
    end
        

end

save overprecfoll_withkeystrokeintersections.mat;

clearvars -except tests *clip *int;
save overprecfoll_formeans.mat;

genocodes = [3 3 3 1 1 1];


% Time Bin 1

over1mn(1,1:length(keys1int)) = mean(cat(1,overrate1_clip{genocodes==1}));
over1mn(2,1:length(keys1int)) = mean(cat(1,overrate1_clip{genocodes==3}));
prec1mn(1,1:length(keys1int)) = mean(cat(1,precrate1_clip{genocodes==1}));
prec1mn(2,1:length(keys1int)) = mean(cat(1,precrate1_clip{genocodes==3}));
foll1mn(1,1:length(keys1int)) = mean(cat(1,follrate1_clip{genocodes==1}));
foll1mn(2,1:length(keys1int)) = mean(cat(1,follrate1_clip{genocodes==3}));

over1sem(1,1:length(keys1int)) = std(cat(1,overrate1_clip{genocodes==1}))/sqrt(length(tests));
over1sem(2,1:length(keys1int)) = std(cat(1,overrate1_clip{genocodes==3}))/sqrt(length(tests));
prec1sem(1,1:length(keys1int)) = std(cat(1,precrate1_clip{genocodes==1}))/sqrt(length(tests));
prec1sem(2,1:length(keys1int)) = std(cat(1,precrate1_clip{genocodes==3}))/sqrt(length(tests));
foll1sem(1,1:length(keys1int)) = std(cat(1,follrate1_clip{genocodes==1}))/sqrt(length(tests));
foll1sem(2,1:length(keys1int)) = std(cat(1,follrate1_clip{genocodes==3}))/sqrt(length(tests));

% Time Bin 2
over2mn(1,1:length(keys2int)) = mean(cat(1,overrate2_clip{genocodes==1}));
over2mn(2,1:length(keys2int)) = mean(cat(1,overrate2_clip{genocodes==3}));
prec2mn(1,1:length(keys2int)) = mean(cat(1,precrate2_clip{genocodes==1}));
prec2mn(2,1:length(keys2int)) = mean(cat(1,precrate2_clip{genocodes==3}));
foll2mn(1,1:length(keys2int)) = mean(cat(1,follrate2_clip{genocodes==1}));
foll2mn(2,1:length(keys2int)) = mean(cat(1,follrate2_clip{genocodes==3}));

over2sem(1,1:length(keys2int)) = std(cat(1,overrate2_clip{genocodes==1}))/sqrt(length(tests));
over2sem(2,1:length(keys2int)) = std(cat(1,overrate2_clip{genocodes==3}))/sqrt(length(tests));
prec2sem(1,1:length(keys2int)) = std(cat(1,precrate2_clip{genocodes==1}))/sqrt(length(tests));
prec2sem(2,1:length(keys2int)) = std(cat(1,precrate2_clip{genocodes==3}))/sqrt(length(tests));
foll2sem(1,1:length(keys2int)) = std(cat(1,follrate2_clip{genocodes==1}))/sqrt(length(tests));
foll2sem(2,1:length(keys2int)) = std(cat(1,follrate2_clip{genocodes==3}))/sqrt(length(tests));

% Time Bin 3
over3mn(1,1:length(keys3int)) = mean(cat(1,overrate3_clip{genocodes==1}));
over3mn(2,1:length(keys3int)) = mean(cat(1,overrate3_clip{genocodes==3}));
prec3mn(1,1:length(keys3int)) = mean(cat(1,precrate3_clip{genocodes==1}));
prec3mn(2,1:length(keys3int)) = mean(cat(1,precrate3_clip{genocodes==3}));
foll3mn(1,1:length(keys3int)) = mean(cat(1,follrate3_clip{genocodes==1}));
foll3mn(2,1:length(keys3int)) = mean(cat(1,follrate3_clip{genocodes==3}));

over3sem(1,1:length(keys3int)) = std(cat(1,overrate3_clip{genocodes==1}))/sqrt(length(tests));
over3sem(2,1:length(keys3int)) = std(cat(1,overrate3_clip{genocodes==3}))/sqrt(length(tests));
prec3sem(1,1:length(keys3int)) = std(cat(1,precrate3_clip{genocodes==1}))/sqrt(length(tests));
prec3sem(2,1:length(keys3int)) = std(cat(1,precrate3_clip{genocodes==3}))/sqrt(length(tests));
foll3sem(1,1:length(keys3int)) = std(cat(1,follrate3_clip{genocodes==1}))/sqrt(length(tests));
foll3sem(2,1:length(keys3int)) = std(cat(1,follrate3_clip{genocodes==3}))/sqrt(length(tests));




save workspace_withmeansandsemsofovprecfoll.mat;
clearvars -except tests *clip *int;


load workspace_withmeansandsemsofovprecfoll.mat;
clearvars -except tests *clip *int;

% QUESTION: Is voc before a chase always higher than during or following?

% Time Bin 1

for i = 1:length(tests) 
    chase(i,1) = precrate1_clip{i}(1); % d
    chase(i,2) = overrate1_clip{i}(1);
    chase(i,3) = follrate1_clip{i}(1);
    selfgroom(i,1) = precrate1_clip{i}(2); % e
    selfgroom(i,2) = overrate1_clip{i}(2);
    selfgroom(i,3) = follrate1_clip{i}(2);
    sniff(i,1) = precrate1_clip{i}(3); % j
    sniff(i,2) = overrate1_clip{i}(3);
    sniff(i,3) = follrate1_clip{i}(3);
    dig(i,1) = precrate1_clip{i}(4); % k
    dig(i,2) = overrate1_clip{i}(4);
    dig(i,3) = follrate1_clip{i}(4);
    explore(i,1) = precrate1_clip{i}(5); % l
    explore(i,2) = overrate1_clip{i}(5);
    explore(i,3) = follrate1_clip{i}(5);
end

graphlabels = {'chase','selfgroom','sniff','dig','explore'};
ticklabels = {'before','during','after'};
graphs = {chase, selfgroom, sniff, dig, explore};

% OK, let's make all the graphs:
for i = 1:length(graphlabels)
figure;
bar(mean(graphs{i}));
hold on;
    for j = 1:length(ticklabels)
    semgraph(j,mean(graphs{i}(:,j)),std(graphs{i}(:,j))./sqrt(length(tests)),'black');
    end
title(graphlabels{i});
ylabel('msec voc/sec behavior');
axis([0 4 0 150]);

end

% Time Bin 2


for i = 1:length(tests) 
    chase(i,1) = precrate2_clip{i}(1); % d
    chase(i,2) = overrate2_clip{i}(1);
    chase(i,3) = follrate2_clip{i}(1);
    selfgroom(i,1) = precrate2_clip{i}(2); % e
    selfgroom(i,2) = overrate2_clip{i}(2);
    selfgroom(i,3) = follrate2_clip{i}(2);
    sniff(i,1) = precrate2_clip{i}(3); % j
    sniff(i,2) = overrate2_clip{i}(3);
    sniff(i,3) = follrate2_clip{i}(3);
    dig(i,1) = precrate2_clip{i}(4); % k
    dig(i,2) = overrate2_clip{i}(4);
    dig(i,3) = follrate2_clip{i}(4);
    explore(i,1) = precrate2_clip{i}(5); % l
    explore(i,2) = overrate2_clip{i}(5);
    explore(i,3) = follrate2_clip{i}(5);
end

graphlabels = {'chase','selfgroom','sniff','dig','explore'};
ticklabels = {'before','during','after'};
graphs = {chase, selfgroom, sniff, dig, explore};

% OK, let's make all the graphs:
for i = 1:length(graphlabels)
figure;
bar(mean(graphs{i}));
hold on;
    for j = 1:length(ticklabels)
    semgraph(j,mean(graphs{i}(:,j)),std(graphs{i}(:,j))./sqrt(length(tests)),'black');
    end
title(graphlabels{i});
ylabel('msec voc/sec behavior');
axis([0 4 0 150]);

end

% Time Bin 3


for i = 1:length(tests) 
    chase(i,1) = precrate3_clip{i}(1); % d
    chase(i,2) = overrate3_clip{i}(1);
    chase(i,3) = follrate3_clip{i}(1);
    selfgroom(i,1) = precrate3_clip{i}(2); % e
    selfgroom(i,2) = overrate3_clip{i}(2);
    selfgroom(i,3) = follrate3_clip{i}(2);
    sniff(i,1) = precrate3_clip{i}(3); % j
    sniff(i,2) = overrate3_clip{i}(3);
    sniff(i,3) = follrate3_clip{i}(3);
    dig(i,1) = precrate3_clip{i}(4); % k
    dig(i,2) = overrate3_clip{i}(4);
    dig(i,3) = follrate3_clip{i}(4);
    explore(i,1) = precrate3_clip{i}(5); % l
    explore(i,2) = overrate3_clip{i}(5);
    explore(i,3) = follrate3_clip{i}(5);
end

graphlabels = {'chase','selfgroom','sniff','dig','explore'};
ticklabels = {'before','during','after'};
graphs = {chase, selfgroom, sniff, dig, explore};

% OK, let's make all the graphs:
for i = 1:length(graphlabels)
figure;
bar(mean(graphs{i}));
hold on;
    for j = 1:length(ticklabels)
    semgraph(j,mean(graphs{i}(:,j)),std(graphs{i}(:,j))./sqrt(length(tests)),'black');
    end
title(graphlabels{i});
ylabel('msec voc/sec behavior');
axis([0 4 0 150]);

end



% Susan wants total vocalization time and number of whistles in each bin so
% let's make a quick spreadsheet for her

% BEGIN SPREADSHEET WRITE

ids = cat(1,tests(:).id);
vtime1 = zeros(1,numel(tests));
vtime2 = zeros(1,numel(tests));
vtime3 = zeros(1,numel(tests));
usvs1 = zeros(1,numel(tests));
usvs2 = zeros(1,numel(tests));
usvs3 = zeros(1,numel(tests));
for i = 1:numel(tests)
    vtime1(i) = tests(i).totvtime1;
    vtime2(i) = tests(i).totvtime2;
    vtime3(i) = tests(i).totvtime3;
    usvs1(i) = size(tests(i).twhis1,2);
    usvs2(i) = size(tests(i).twhis2,2);
    usvs3(i) = size(tests(i).twhis3,2);
end

vtimes = [ids,vtime1', vtime2', vtime3'];
usvs = [ids, usvs1', usvs2', usvs3'];

xlswrite('krystalc6_usvs.xls',vtimes,'totvtime','A1');
xlswrite('krystalc6_usvs.xls',usvs,'totUSVs','A1');
clear ids vtime* usv*;




