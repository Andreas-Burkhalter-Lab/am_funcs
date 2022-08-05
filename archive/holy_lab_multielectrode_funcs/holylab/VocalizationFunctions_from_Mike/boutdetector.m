function [allBouts] = boutdetector(data,cutoff)

% This code assumes that directory has been trunced to include only
% sngs from animals with >= 10 whistles per sng


a = (size(data,1)-1);
for i=1:a
twhises{i} = data{i+1,find(strcmp('twhis',data(1,:)))};
end


n = size(twhises,2); % number of animals to loop over
allBouts = cell(1,n);
    for animal = 1:n
        twhis = twhises{animal};
        if size(twhis,2) >= 2
            
        whis1 = twhis;
        whis2 = twhis;
        
        whis1 = whis1(2,:); % end times of the whistles
        whis2 = whis2(1,:); % start times of the whistles
        
        whis1 = whis1(1:length(whis1)-1); %end times of whistle x
        whis2 = whis2(2:length(whis2)); %start times of whistle x+1
        
        % make a vector of pauses
           
        pauses = whis2 - whis1; % pauses(1) means twhis(1,2) - twhis(2,1)
        
        boundaries = find(pauses >= cutoff); % these ids correspond to the ids of whis1
        
        starts = [1 (boundaries+1)]; % starts go from id 1 in twhis through each id in boundaries +1
        ends = [boundaries size(twhis,2)]; %ends go from id 1 in boundaries through last id in boundaries
        
                
        boutstarttimes = twhis(1,starts);
        boutendtimes = twhis(2,ends);
        whiscounts = ends - starts;
        whiscounts = whiscounts + 1; % if end is 4 and start is 1, then it will say there are 3 instead of 4.
                                    % this leaves zeros too, which is
                                    % evidence of this. So add 1 to all the
                                    % counts
        
        boutinfo = [boutstarttimes; boutendtimes;whiscounts];
        allBouts{animal} = boutinfo;
        
          
        clear boutstarttimes boutendtimes whiscounts pauses boundaries starts ends whis1 whis2 twhis;
        elseif size(twhis,2) == 1
        allBouts{animal} = [twhis(1,1);twhis(2,1);1];
        else
        allBouts{animal} = [];
        end

    end








end
