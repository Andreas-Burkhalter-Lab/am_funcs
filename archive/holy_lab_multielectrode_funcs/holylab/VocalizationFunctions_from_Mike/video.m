prefix = input('Prefix? ','s');
suffix = input('Suffix? ','s');
wildcard = ['*',suffix];
a = dir(wildcard);
csvfile_id = zeros(1,numel(a));
% Now let's make a vector of the ids 
%
for i = 1:numel(a);
    n=a(i).name;
    n=strrep(n,suffix,'');
    n=strrep(n,prefix,'');
    % now it's of the form ID
    n=cat(2,n);
    csvfile_id(i) = str2num(n);
end
  
clearvars -except a csvfile_id;
%%
%I.  READ THE CSVS.
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

clearvars -except a csv*

%%
% ADD THE CSVS TO A TABLE
videodata = {'id','whisn','twhis','whispf','whistype','csvs','keynames','keyscores','whisdt','sum(whisdt)','keydt','sum(keydt)','oversum','overidx','overrate','precsum','precidx','precrate','follsum','follidx','follrate','keyintersect','intersectrate'};
for i = 1:length(csvfile_id)
    videodata{i+1,1} = csvfile_id(i);
    videodata{i+1,6} = csvs{i};
end

for i = 1:(size(videodata,1)-1)
        keys = videodata{i+1,6}(:,1);
        startstop = videodata{i+1,6}(:,2);
        times = videodata{i+1,6}(:,3);
        keys = keys';
        startstop = startstop';
        times = times';
            for j = 1:length(times)
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
                % Clip extra starts and stops
                x = [x;0]; %add a spacer
                x1 = x(1:end-1);
                x2 = x(2:end);
                noextra = find(x1 ~= x2);
                x = x(noextra);
            
            t = timesdouble(k);
                % Clip extra start and stop times
                t = [t 0]; % add a spacer
                t = t(noextra);
                
            y = [t(find(x==1));t(find(x == 2))];
            %y(3,1:size(y,2)) = j;
            % There is still one more thing to clean up. Apparently
            % sometimes start time = stop time. I don't know how this 
            % happens but assume it's a fuckup.
            [~,c] = find(y(1,:) ~= y(2,:)); % Makes sure all starts are 
            % different from stops
            y = y(:,c);
            keyscores{j} = y;
            clear k x* t y noextra c;
        end  
         % Now, double check that every keyname has a non-empty bout.
         % Because of the cleanups (starts without stops or vice versa, or
         % starts and stops that are the same), there may end up being keys
         % that have no score.
         for j = 1:length(keyscores)
             if ~isempty(keyscores{j})
                 keepers(j) = j;
             else
                 keepers(j) = NaN;
             end
         end
         keepers = keepers(~isnan(keepers));
         keyscores = keyscores(keepers);
         keynames = keynames(keepers);
         videodata{i+1,7} = keynames;
         videodata{i+1,8} = keyscores;
         clear key* start* times* keepers;
end
%%
% 1:        twhis start or end within keyscore start&end OVERLAP
% 2:        twhis end within 5 sec of keyscore start     PRECEDE
% 3:        twhis start within 5 sec of keyscore end     FOLLOW
% 0:        none of the above

% Calculate whistle duration times and total vocalization time
%vocfile = input('Name of file with vocalization data?','s');
%load(['../Vars/',vocfile]);
clearvars -except videodata data;
for i = 1:(size(videodata,1)-1)
    videodata{i+1,2} = data{i+1,6};
    videodata{i+1,3} = data{i+1,7};
    videodata{i+1,4} = data{i+1,28};
    videodata{i+1,5} = data{i+1,32};
end
for i = 1:(size(videodata,1)-1)
  if ~isempty(videodata{i+1,2})
    videodata{i+1,9} = videodata{i+1,3}(2,:) - videodata{i+1,3}(1,:);
    videodata{i+1,10} = sum(videodata{i+1,9});
  end
end


% Calculate behavior duration times

for i = 1:(size(videodata,1)-1)
   if ~isempty(videodata{i+1,7})
       tkey = videodata{i+1,8};
       keydt = cell(1,length(tkey));
       keydtsum = zeros(1,length(tkey));       
       for j = 1:length(tkey)
        tkey{j} = tkey{j}(1:2,:);
        keydt{j} = tkey{j}(2,:) - tkey{j}(1,:);
        keydtsum(j) = sum(keydt{j});
       end
       videodata{i+1,11} = keydt;
       videodata{i+1,12} = keydtsum; 
   end
end
%%


%Score Vocalization/Keystroke Overlap
for i = 1:(size(videodata,1)-1)
    if ~isempty(videodata{i+1,5}) % Make sure there are scores
    
        x = videodata{i+1,8}; %tkey
        twhis = videodata{i+1,3}; %twhis
        dt = videodata{i+1,9}; %whisdt
        
        for key = 1:length(x);
            y = x{key};
            
            for bout = 1:size(y,2) % loop over each bout for keyth key
            
                bstart = y(1,bout);
                bend = y(2,bout);
                
                
                overlaps = find((twhis(1,:)>= bstart & twhis(2,:) <= bend) | ...
                                (twhis(1,:)<= bstart & twhis(2,:) >= bstart)| ...
                                (twhis(1,:)<= bend & twhis(2,:) >= bend));
                [~,overindx] = find((twhis(1,:)>= bstart & twhis(2,:) <= bend) | ...
                                (twhis(1,:)<= bstart & twhis(2,:) >= bstart)| ...
                                (twhis(1,:)<= bend & twhis(2,:) >= bend));
                            
                precedes = find((twhis(1,:)<= bstart & twhis(2,:) <= bstart) & ...
                                (abs(twhis(2,:) - bstart) <= 5));
               [~,precindx] = find((twhis(1,:)<= bstart & twhis(2,:) <= bstart) & ...
                                (abs(twhis(2,:) - bstart) <= 5));
                            
                follows = find((twhis(1,:) >= bend) & (twhis(2,:) >= bend) & ...
                                (abs(twhis(1,:) - bstart) <= 1));
               [~,follindx] = find((twhis(1,:) >= bend) & (twhis(2,:) >= bend) & ...
                                (abs(twhis(1,:) - bstart) <= 5));            
                            
                oversum{key}(bout) = sum(dt(overlaps));
                precedesum{key}(bout) = sum(dt(precedes));
                followsum{key}(bout) = sum(dt(follows));
                overix{key}{bout} = overindx;
                precix{key}{bout} = precindx;
                follix{key}{bout} = follindx;
            
                clear bstart bend overlaps precedes follows;
                clear *indx;
            end
            clear y;
        end
        
        clear x twhis dt;
        
        overs{i} = oversum;
        precs{i} = precedesum;
        folls{i} = followsum;
        ox{i} = overix;
        px{i} = precix;
        fx{i} = follix;
        
        clear oversum precedesum followsum overix precix follix;
    end
end

% Calculate Rate msec Voc/ sec Behavior
for i = 1:numel(overs)
    for key = 1:numel(overs{i})
        keydt = videodata{i+1,12}(key);
        overrate{i}(key) = (sum(overs{i}{key})*1000)/keydt;
        precrate{i}(key) = (sum(precs{i}{key})*1000)/keydt;
        follrate{i}(key) = (sum(folls{i}{key})*1000)/keydt;
    end
end

clearvars -except videodata over* prec* foll* ox px fx;

for i = 1:(size(videodata,1)-1)
    videodata{i+1,13} = overs{i};
    videodata{i+1,14} = ox{i};
    videodata{i+1,15} = overrate{i};
    videodata{i+1,16} = precs{i};
    videodata{i+1,17} = px{i};
    videodata{i+1,18} = precrate{i};
    videodata{i+1,19} = folls{i};
    videodata{i+1,20} = fx{i};
    videodata{i+1,21} = follrate{i};
end

clearvars -except videodata;

% Find the instersection of keys scored in all files in each bin.
%%

for i = 1:(size(videodata,1)-1)
    keys{i} = cat(2,videodata{i+1,7}{:});
end

keysint = keys{1};
for i = 2:length(keys)
    keysint = intersect(keysint,keys{i});
end
save('../Vars/keysint.mat','keysint');
for i = 1:(size(videodata,1)-1)
    videodata{i+1,22} = keysint;
end

for i = 1:(size(videodata,1)-1)
    overrate = videodata{i+1,15};
    precrate = videodata{i+1,18};
    follrate = videodata{i+1,21};
    intrate = cell(1,length(keysint));
    for key = 1:length(keysint)
        idx = find(keys{i} == keysint(key));
        intrate{key} = [precrate(idx), overrate(idx), follrate(idx)];
        clear idx;
    end
    videodata{i+1,23} = intrate;    
    clear overrate precrate follrate intrate;
end

%% OVERLAY PLOT CODE
clearvars -except videodata;
% Overlay behaviors onto the SNG plot
load('../Vars/keymeanings.mat');
% Plot the sonogram:
    sngfile = input('What id# ?');
% Grab the keys for this file and make list of meanings
    keys = videodata{find(cat(2,videodata{2:end,1})==(sngfile))+1,7};
    list = cell(1,length(keys));
    for i = 1:length(keys)
        list{i} = keymeanings(find(strcmp(keys{i},cat(2,keymeanings(:,2)))),1);
    end
% Grab start and stop times
    tkey = videodata{find(cat(2,videodata{2:end,1})==(sngfile))+1,8};
% Y-spacers for the graph
    yspacer = zeros(1,length(keys));
    yspacer(1) = 50;
    for i = 2:length(yspacer)
        yspacer(i)=yspacer(i-1)-2;
    end
% Color codes
    x = colormap('HSV');
    close;
    n = randn(1,size(x,1));
    x = [x,n'];
    x = sortrows(x,4);
    x = x(:,1:3);
    clear n;
    x = x(1:length(yspacer),:);
% Plot tkeys for each color and spacer

for i = 1:length(list)
    table{i} = strcat(num2str(i),'. ',list{i});
end
table = cat(1,table{:});
disp(table);
howmany = input('Plot which keys?');
whatkind = input('Plot (1) Sonogram or (2) twhis/meanfreq? ');
    if whatkind == 1
     spsngplot(strcat('cMG2013_',num2str(sngfile,'%03d'),'.WAV.SNG'));
    elseif whatkind ==2
     twhis = cat(2,videodata{find(cat(2,videodata{2:end,1})==(sngfile))+1,3});
     pf = cat(2,videodata{find(cat(2,videodata{2:end,1})==(sngfile))+1,4});
     mf = zeros(1,length(pf));
      for i = 1:length(pf)    
        mf(1:2,i) = mean(pf{i})/1000;
      end
     plot(twhis,mf,'red','LineWidth',2);
    end
xlabel('time ( sec ) ');
ylabel('frequency ( kHz ) ');
titletext = {(num2str(sngfile,'%03d'))};
    for i = 1:length(howmany)
        titletext(i+1) = strcat(num2str(i),'.',list{howmany(i)});
    end
titletext = titletext';    
title(titletext);
hold on;
plot(tkey{howmany(1)},[yspacer(1) yspacer(1)],'Color',x(1,1:3),'LineWidth',3);
for i = 2:length(howmany)
    plot(tkey{howmany(i)},[yspacer(i) yspacer(i)],'Color',x(i,1:3),'LineWidth',3);
end

%set(gca,'YTick',[0,20,30:2:50,60:20:120]);

% END OVERLAY PLOT


%%
% Probability that whis is in any of the shared keystrokes

clearvars -except videodata;

for i = 1:(size(videodata,1)-1)
overidx = videodata{i+1,14};
whisn = videodata{i+1,2};
keyname = videodata{i+1,7};
keysect = videodata{i+1,22};
keydt = videodata{i+1,12};

    for j = 1:length(keysect)
    over = overidx{find(strcmp(keysect(j),keyname))};
    lengths = zeros(1,length(over));
    dt = keydt(find(strcmp(keysect(j),keyname)));
       for count = 1:length(over)
           lengths(count) = length(over{count});
       end
    videodata{i+1,23+j}  = (sum(lengths)/whisn)/(dt/60);
    clear idx over lengths dt;
    end
clear overidx whisn keyname keysect keydt;
end
%% (whisn,meanfreq,f(ss),f(jumps) associated with each behavior
%clearvars -except videodata;
load('../Vars/keymeanings.mat');
adder = 30;%size(videodata,2);
%firstcol = size(videodata,2) + 1;
%lastcol = (size(keymeanings,1)-1) + firstcol;

for i = 1:(size(videodata,1)-1)
    overidx = videodata{i+1,16};%14};
    whisn = videodata{i+1,4};%2};
    whispf = videodata{i+1,6};%4};
    whistype = videodata{i+1,7};%5};
    keyname = videodata{i+1,9};%7};
    sumkeydt = videodata{i+1,14};
    
        for j = 2:size(keymeanings,1)
        keyidx = find(strcmp(keymeanings(j,2),keyname));
            if isempty(keyidx)
                videodata{i+1,adder+j-1} = [NaN,NaN,NaN,NaN];
            else
                ox = overidx{keyidx}; 
                keydt = sumkeydt(keyidx);
                ox = cat(2,ox{:});
                vec(1) = length(ox)/keydt; % # whistles/sec behavior
                
                pfs = whispf(ox);
                vec(2) = (mean(cat(2,pfs{:})))/1000;
                types = whistype(ox);
                ssidx = find(strcmp('s',types));
                jumpidx = find(~strcmp('s',types));
                vec(3) = length(ssidx)/length(ox);
                vec(4) = length(jumpidx)/length(ox);
                
                videodata{i+1,adder+j-1} = vec;
                
                clear ox keydt vec pfs types ssidx jumpidx;
            end
        end
    
end
% Column headers
% for j = 2:size(keymeanings,1)
%     videodata{1,adder+j-1} = keymeanings{j,1};
% end


%%
% Number & Duration of Bouts of Mounting
clearvars -except videodata;
load('../Vars/keymeanings.mat');
[r,~] = find(strcmp('Mount',keymeanings));
adder = size(videodata,2);
for i = 1:(size(videodata,1)-1)
    keyname = videodata{i+1,9};
    keyscores = videodata{i+1,10};
    keydt = videodata{i+1,13};
    
    mountidx = find(strcmp(keymeanings(r,2),keyname));
    
    if ~isempty(mountidx)
        keyscores = keyscores{mountidx};
        keydt = keydt{mountidx};
        videodata{i+1,adder+1} = size(keyscores,2);
        videodata{i+1,adder+2} = mean(keydt);    
    else
        videodata{i+1,adder+1} = 0;
        videodata{i+1,adder+2} = 0;  
    end
end





%%
% EDIT EVERYTHING BELOW

for i = 1:18
asniff(i) = videodata{i+1,29}(4);
tsniff(i) = videodata{i+1,30}(4);
hsniff(i) = videodata{i+1,31}(4);
sgroom(i) = videodata{i+1,33}(4);
follow(i) = videodata{i+1,36}(4);
explore(i) = videodata{i+1,38}(4);
mount(i) = videodata{i+1,40}(4);
end
% Look at just anogenital sniffing over bouts.
clearvars -except videodata;
figure;
%hold on;
keyinput = input('What key do you want? ','s');
for i = 1:(size(videodata,1)-1)
    keys = videodata{i+1,5};
    idx = find(strcmp(keys,keyinput)==1);
    over = videodata{i+1,11};
    geno = videodata{i+1,24};
    forplot = over{idx};
    if geno == 1
        plot(forplot,'blue');
    elseif geno ==3
        plot(forplot,'red');
    end
    pause;
    disp('Hit enter for the next one');
    clear keys idx over geno forplot;
end

% Make a table who is mxn animals x max # of bouts of keystroke

clearvars -except videodata;

keyinput = input('What key do you want? ','s');
for i = 1:(size(videodata,1)-1)
    keys = videodata{i+1,5};
    idx = find(strcmp(keys,keyinput)==1);
    over = videodata{i+1,11};
    boutdata{i} = over{idx};
    nbout(i) = size(over{idx},2);
end

table = zeros((size(videodata,1)-1),max(nbout));
for i = 1:length(boutdata)
    z = boutdata{i};
    for j = 1:size(z,2)
        table(i,j) = z(j);
    end
    clear z;
end

plot(mean(table));
xlabel('bout of behavior');
ylabel('amount of vocalization (sec)');
axis([0 60 0 3]);

%% Percent of Bout With Vocalization.

% Requires the column keydt, that has the duration time for each bout of
% each scored behavior, and the column oversum, that has the duration time
% for vocalization during each bout of each behavior.


for i = 2:size(videodata,1)
    keydt = videodata{i,find(strcmp('keydt',cat(2,videodata(1,1:end))))};
    oversum = videodata{i,find(strcmp('oversum',cat(2,videodata(1,1:end))))};
    coverage = cell(1,length(keydt));
    for j = 1:length(keydt)
        fraction = oversum{j}./keydt{j};
        % Clips Infs and NaN -
        % Infs occur because occasionally a key was hit and immediately
        % released, resulting in an error where a bout has dt = 0
        % This didn't effect previous analyses because vocalization cannot
        % co-occur with 0 so it was ignored. Now that I'm doing the
        % opposite sort of operation (bouts occuring with vocalization), it
        % shows up. In the data set there are very few of these, but they
        % do occur.
        fraction = fraction(find(~isnan(fraction)));
        fraction = fraction(find(~isinf(fraction)));
        % Clip all fraction greater than or equal to 1. This occurs because
        % the oversum includes whistles that overlap the boundaries. If the
        % bout itself is shorter than the whistle that overlaps it, the
        % fraction is greater than 1. But these are probably errors because
        % these would be bouts of behaviors ~ 100 msec long which can
        % probably not be reliably scored.
        fraction = fraction(find(fraction < 1));
        videodata{i,50}(j) = length((find(fraction>=0.1)))/length(fraction);
        coverage{j} = fraction;
        clear fraction;
    end
    videodata{i,49} = coverage;
    clear coverage;
    
end

%% Distribution of bouts with and without vocalization
clearvars -except allvideo
load('../females/Vars/keymeanings.mat');
for i=2:size(keymeanings,1)
disptext{i} = ...
    [num2str(i-1),' . (',keymeanings{i,2},') ',keymeanings{i,1}];
end
disptext = disptext';
disp(disptext(2:end));
key = input('Which key? Enter the key : ','s');
[keyrow,~] = find(strcmp(key,keymeanings));
keydescription = keymeanings{keyrow,1};
% Now, based on the key, plot the distribution:

for i = 2:size(allvideo,1)
    keys = allvideo{i,find(strcmp('keynames',cat(2,allvideo(1,1:end))))};
    keyidx = find(strcmp(key,keys));
    if ~isempty(keyidx)
    coverage = allvideo{i,find(strcmp('coverage',cat(2,allvideo(1,1:end))))};
    x{i-1} = coverage{keyidx};
    end
end

x = cat(2,x{:});

bins = 0:0.05:1;
h = hist(x,bins);
relh = (h./sum(h))*100;
figure;
bar(bins,relh);
title({['Vocalization during bouts of ',keydescription];...
      ['Number of bouts scored = ',num2str(sum(h))];...
      ['Bin width = ', num2str(bins(2)-bins(1))];...
      });
      
 xlabel('Fraction of bout covered by vocalization');
 ylabel('% of all bouts');
 axis([0 1 0 100]);


%% Distribution of fraction of bouts with greater than 10% vocalization

clearvars -except allvideo
load('../females/Vars/keymeanings.mat');
for i=2:size(keymeanings,1)
disptext{i} = ...
    [num2str(i-1),' . (',keymeanings{i,2},') ',keymeanings{i,1}];
end
disptext = disptext';
disp(disptext(2:end));
key = input('Which key? Enter the key : ','s');
[keyrow,~] = find(strcmp(key,keymeanings));
keydescription = keymeanings{keyrow,1};
% Now, based on the key, plot the distribution:

for i = 2:size(allvideo,1)
    keys = allvideo{i,find(strcmp('keynames',cat(2,allvideo(1,1:end))))};
    keyidx = find(strcmp(key,keys));
    if ~isempty(keyidx)
    withvoc = allvideo{i,find(strcmp('withvoc',cat(2,allvideo(1,1:end))))};
    x(i-1) = withvoc(keyidx);
    end
end

x = cat(2,x{:});

bins = 0:0.05:1;
h = hist(x,bins);
relh = (h./sum(h))*100;
figure;
bar(bins,relh);
title({['Vocalization during bouts of ',keydescription];...
      ['Number of bouts scored = ',num2str(sum(h))];...
      ['Bin width = ', num2str(bins(2)-bins(1))];...
      });
      
 xlabel('Fraction of bout covered by vocalization');
 ylabel('% of all bouts');
 axis([0 1 0 100]);
