%% Make a list of IDs
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

%% Read CSVs

for i = 1:numel(a)
    x = readacsv(a(i).name);
    csvs{i} = x;
end

clearvars -except csv*;

%% Load CSVs and IDs into data table
videodata = {'id','genotype','sex','csv'};

for i = 1:length(csvfile_id)
    videodata{i+1,find(strcmp(cat(2,videodata(1,:)),'id'))} = csvfile_id(i);
    videodata{i+1,find(strcmp(cat(2,videodata(1,:)),'csv'))} = csvs{i};
end

for i = 1:length(behaviornames)
    videodata{1,find(strcmp(cat(2,videodata(1,:)),'csv'))+i} = behaviornames{i};
end

%% Split CSVs
clearvars -except videodata;
for i = 1:size(videodata,1)-1
    csv = videodata{i+1,4};
    [~,subject] = find(strcmp(' Subject ',csv) | strcmp('Subject',csv));
    [~,state] = find(strcmp(' Event_Type ',csv) | strcmp('Event_Type',csv));
    [~,behavior] = find(strcmp(' Behavior ',csv)| strcmp('Behavior',csv));
    [~,timestamp] = find(strcmp(' Time_Relative_sf ',csv)| strcmp('Time_Relative_sf',csv));
    
    subs{i} = csv(2:end,subject);
    states{i} = csv(2:end,state);
    behaviors{i} = csv(2:end,behavior);
    timestamps{i} = csv(2:end,timestamp);
    
    for j = 1:length(subs{i})
        subs{i}{j} = strrep(subs{i}{j},' ','');
        states{i}{j} = strrep(states{i}{j},' ','');
        behaviors{i}{j} = strrep(behaviors{i}{j},' ','');
        timestamps{i}{j} = strrep(timestamps{i}{j},' ','');
        times{i}(j) = str2num(timestamps{i}{j});
    end
end
clearvars -except videodata times subs behaviors states; 

for i = 1:size(videodata,1)-1
   
        b = behaviors{i};
        s = states{i};
        t = times{i};
        c = subs{i};
    
        [testid,~] = find(strcmp(c,'Testanimal'));
        [stimid,~] = find(strcmp(c,'Stimulusanimal'));
        
        b_stim = b(stimid);
        b_test = b(testid);
        s_stim = s(stimid);
        s_test = s(testid);
        t_stim = t(stimid);
        t_test = t(testid);
        
    for j = 5:15
        [name_test,~] = find(strcmp(b_test,videodata{1,j}));
        if ~isempty(name_test)
        startstop = s_test(name_test);
        u = t_test(name_test);
        start = find(strcmp(startstop,'Statestart'));
        stop = find(strcmp(startstop,'Statestop'));
        bouts{1}(1,1:length(start)) = u(start);
        bouts{1}(2,1:length(start)) = u(stop);
            if bouts{1}(1,1) == 0 && bouts{1}(2,1) ==0
                bouts{1} = bouts{1}(1:2,2:end);
            end
        bouts{2} = bouts{1}(2,:) - bouts{1}(1,:);
        bouts{3} = size(bouts{1},2);
        bouts{4} = mean(bouts{2});
        clear name_test;
        clear startstop;
        clear start stop;
        else
            bouts{1} = [0,0];
            bouts{2} = 0;
            bouts{3} = 0;
            bouts{4} = 0;
        end
        videodata{i+1,j} = bouts;
        clear bouts;
    end
    
    for j = 16:26
        [name_stim,~] = find(strcmp(b_stim,videodata{1,j}));
        if ~isempty(name_stim)
        startstop = s_stim(name_stim);
        u = t_stim(name_stim);
        start = find(strcmp(startstop,'Statestart'));
        stop = find(strcmp(startstop,'Statestop'));
        bouts{1}(1,1:length(start)) = u(start);
        bouts{1}(2,1:length(start)) = u(stop);
            if bouts{1}(1,1) == 0 && bouts{1}(2,1) ==0
                bouts{1} = bouts{1}(1:2,2:end);
            end
        bouts{2} = bouts{1}(2,:) - bouts{1}(1,:);
        bouts{3} = size(bouts{1},2);
        bouts{4} = mean(bouts{2});
        clear name_test;
        clear startstop;
        clear start stop;
        else
            bouts{1} = [0;0];
            bouts{2} = 0;
            bouts{3} = 0;
            bouts{4} = 0;
        end
        videodata{i+1,j} = bouts;
        clear bouts;
    end
end
clearvars -except videodata;
%% Behavior table for stats

table = videodata(1:end,1:3);
tests = cat(2,videodata(1,5:15));
sniffs = [find(strcmp(tests,'Anogenitalsniffing')),...
                find(strcmp(tests,'Headsniffing')),...
                find(strcmp(tests,'Trunksniffing'))];        
grooms = [find(strcmp(tests,'Self-grooming')),...
                find(strcmp(tests,'Allo-grooming'))];
moves = [find(strcmp(tests,'Cageexploration')),...
                find(strcmp(tests,'Chasing')),...
                find(strcmp(tests,'Fleeing'))];
for i = 1:size(videodata,1)-1

    testsniff = videodata(i+1,4+sniffs);
    stimsniff = videodata(i+1,4+sniffs+11);
    testgroom = videodata(i+1,4+grooms);
    stimgroom = videodata(i+1,4+grooms+11);
    testmove = videodata(i+1,4+moves);
    stimmove = videodata(i+1,4+moves+11);
    
    for j = 1:length(testsniff)
        table{i+1,4}(j) = testsniff{j}{5};
        table{i+1,5}(j) = stimsniff{j}{5};
    end
    for j = 1:length(testgroom)
        table{i+1,6}(j) = testgroom{j}{5};
        table{i+1,7}(j) = stimgroom{j}{5};
    end
    for j = 1:length(testmove)
        table{i+1,8}(j) = testmove{j}{5};
        table{i+1,9}(j) = stimmove{j}{5};
    end
end
         
%% Move tot voc and genotype based on id

for i = 1:size(table,1)-1
id = table{i+1,1};
k = find(ids ==id);
table{i+1,2} = geno(k);
table{i+1,10} = x(k);
end
%% Overlap calculations

for i = 1:17
    twhis = foroverlap{i+1,2};
    dur = twhis(2,:) - twhis(1,:);
    
    for j = 4:19
        bouttimes = foroverlap{i+1,j}{1};
        boutdur = foroverlap{i+1,j}{2};
        if boutdur ~= 0
        k = [];
        for t = 1:size(bouttimes,2)
        [~,col] = find(twhis(1,:) >= bouttimes(1,t) & twhis(2,:) <= bouttimes(2,t));
        k = [k col];
        clear col;
        end
        if ~isempty(k)
            v = sum(dur(k));
            v = v*1000;
            d = v/sum(boutdur);
        else
            d = 0;
        end
        else
            d = NaN;
        end
        foroverlap{i+1,j}{5} = d;
    end
end

%% Behavior table for stats of vocs

table = videodata(1:end,1:3);
tests = cat(2,videodata(1,5:15));
sniffs = [find(strcmp(tests,'Anogenitalsniffing')),...
                find(strcmp(tests,'Headsniffing')),...
                find(strcmp(tests,'Trunksniffing'))];        
grooms = [find(strcmp(tests,'Self-grooming')),...
                find(strcmp(tests,'Allo-grooming'))];
moves = [find(strcmp(tests,'Cageexploration')),...
                find(strcmp(tests,'Chasing')),...
                find(strcmp(tests,'Fleeing'))];
for i = 1:size(videodata,1)-1

    testsniff = videodata(i+1,4+sniffs);
    stimsniff = videodata(i+1,4+sniffs+11);
    testgroom = videodata(i+1,4+grooms);
    stimgroom = videodata(i+1,4+grooms+11);
    testmove = videodata(i+1,4+moves);
    stimmove = videodata(i+1,4+moves+11);
    
    for j = 1:length(testsniff)
        table{i+1,4}(j) = testsniff{j}{5};
        table{i+1,5}(j) = stimsniff{j}{5};
    end
    for j = 1:length(testgroom)
        table{i+1,6}(j) = testgroom{j}{5};
        table{i+1,7}(j) = stimgroom{j}{5};
    end
    for j = 1:length(testmove)
        table{i+1,8}(j) = testmove{j}{5};
        table{i+1,9}(j) = stimmove{j}{5};
    end
end
         
%% Move tot voc and genotype based on id

for i = 1:size(table,1)-1
id = table{i+1,1};
k = find(ids ==id);
table{i+1,2} = geno(k);
table{i+1,10} = x(k);
end
