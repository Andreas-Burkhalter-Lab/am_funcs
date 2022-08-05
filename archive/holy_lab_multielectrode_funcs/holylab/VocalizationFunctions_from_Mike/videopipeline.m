%% Read csvs

a = dir('*.csv');
csvs = cell(1,numel(a));

for i = 1:numel(a)
    y = readacsv(a(i).name);
    for j = 1:size(y,1)
        if ~isnumeric(y{j,3})
            y{j,3} = str2num(y{j,3});
        end
    end
    csvs{i} = y;
    clear y;
end

clearvars -except csvs;

%% Find key bout times

if ~exist('keys','var')
    load('../Vars/keys.mat')
end

for csv = 1:length(csvs)
x = csvs{csv};
for i = 1:size(keys,1)

    keyidx = find(strcmp(keys{i,1},cat(2,x(1:end,1))));
    
        if ~isempty(keyidx)
            startstop = grp2idx(cat(2,x(keyidx,2)))';
            times = cat(2,x{keyidx,3});


    
    % Cleanups:
    % 1. Remove any excess starts that don't have stops
            if startstop(end) == 1
                startstop = startstop(1:end-1);
                times = times(1:end-1);
            end
            startstop = [startstop,0] ; % adds spacer
            noextra = find(startstop(1:end-1) ~= startstop(2:end));
            startstop = startstop(noextra);
    
    % 2. Remove any bouts < 1 second in duration
            times = times(noextra);
            keytime(1,1:length(times)/2) = times(startstop ==1);
            keytime(2,1:length(times)/2) = times(startstop ==2);
            keytimedur = keytime(2,:) - keytime(1,:);
        
            durid = find(keytimedur >= 1);
            keytime = keytime(1:2,durid);
            
     % 3. Merge bouts with <= 0.5 sec in pause
            if size(keytime,2) > 1
            pauses = keytime(1,2:end) - keytime(2,1:end-1);
            pauseid = find(pauses > 0.5);
            toprow = [1,pauseid+1];
            bottomrow = [toprow(2:end)-1,toprow(end)];
            keytime2(1,:) = keytime(1,toprow);
            keytime2(2,:) = keytime(2,bottomrow);
            else
            keytime2 = keytime;
            end
            
            
            if ~isempty(keytime2)            
            keybouts{i} = keytime2;
            keydurs{i} = keytime2(2,1:end) - keytime2(1,1:end);
            keydursum(i) = sum(keydurs{i});
            keyboutn(i) = length(keybouts{i});
            else
            keybouts{i} = [];
            keydurs{i} = [];
            keydursum(i) = NaN;
            keyboutn(i) = 0;
            end
            clear keyidx startstop times keytime keytimedur durid noextra pauses pauseid toprow bottomrow keytime2;
        else
            keybouts{i} = [];
            keydurs{i} = [];
            keydursum(i) = NaN;
            keyboutn(i) = 0;
        end
end

xbouts{csv} = keybouts;
xdurs{csv} = keydurs;
xdurmsum{csv} = keydursum;
xboutn{csv} = keyboutn;
clear keybouts;
clear keydurs;
clear keydursum;
clear keyboutn;
clear x;


end

%% Make data table

videodata = {'fileid','experiment','animalid','geno','sex','twhis','whisdt','whispf','whistype'};
lastcol = size(videodata,2);

nextcol = lastcol + 1;
for i = 1:size(keys,1)
    videodata{1,nextcol} = keys{i,2};
    videodata{1,nextcol+size(keys,1)} = ['dur ',keys{i,2}];
    videodata{1,nextcol+2*size(keys,1)} = ['total duration ',keys{i,2}];
    videodata{1,nextcol+3*size(keys,1)} = ['# of bouts ' ,keys{i,2}];
    for j = 1:length(xbouts)
        videodata{j+1,nextcol} = xbouts{j}{i};
        videodata{j+1,nextcol+size(keys,1)} = xdurs{j}{i};
        videodata{j+1,nextcol+2*size(keys,1)} = xdurmsum{j}(i);
        videodata{j+1,nextcol+3*size(keys,1)} = xboutn{j}(i);
    end
    nextcol = nextcol+1;
end
clear nextcol lastcol;
