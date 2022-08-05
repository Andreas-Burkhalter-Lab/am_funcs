% This will make a fake sonogram depending on 'subtype' using the Scattoni
% classification scheme.

% For testing round, just load any whistle snip into a variable sng. This
% is just the sparse sonogram for that whistle.

% First, some measurements:

% For 3 minutes of song, this is 3*60 = 180 sec *1000 = 180,000 ms.
% So, we need to make a matrix which is 257 x 180,000

% From my WT brj56 data, the average intrabout pause was 145.8±2.9
% We have a discreet number of columns, so:

% To get a good range, let's take 10 SEMs 

minIntrapause = 10;
maxIntrapause = round(145.8+10*2.9);

% Furthermore, the average interbout pause was 5095.5±1847.4

minInterpause = 260;
maxInterpause = round(5095.5 + 5*1847.4);

% Now, the average number of bouts was 58.7±7.5

minBoutn = 1;
maxBoutn = round(58.7 + 10*7.5);

% The average number of syllables per bout was 2.8±0.2

minBoutlength = 1;
maxBoutlength = round(2.8 + 10*0.2);

% now take a measurement of the size of our example whistle

ncols = size(sng,2);

% Now, we need to decide on a random number of a bunch of things:


% 1) # of bouts & 2) Pauses between them & 3) Boutlengths

Boutn = randi([minBoutn maxBoutn],1);
InterPauses = randi([minInterpause maxInterpause],1,Boutn -1 );
Boutlengths = randi([minBoutlength maxBoutlength],1,Boutn);

% The second column of Boutlengths will now contain the number of pauses
% each bout must possess:

Boutlengths(2,:) = Boutlengths(1,:) - 1 ;

% Great, now construct each Bout in a cell array
% each Bouts{i} = (ncols_whistle, ncols_whistle, ncols_whistle;
%                  (ncols_pause, ncols_pause, 0);


for i = 1:Boutn
    n = Boutlengths(1,i);
 
        for j = 1:n
        z{1,j} = sng;
        z{2,j} = sparse(zeros(257,randi([minIntrapause maxIntrapause],1,1)));
        end
        
    Bouts{i} = z;
    clear z n;
end

% Now we want to sequence these things into one sequence for each bout

for i = 1:size(Bouts,2)
    z = Bouts{i};
    n = numel(z);
    y = {};
    
    for j = 1:n-1
        y{j} = z{j};
    end
    
    Sequenced{i} = y;
    clear z n y;
    
end

% Great, now everything is sequenced, we just need to overlay the pauses

for i = 1:size(Sequenced,2)-1
    Sequenced{2,i} = {sparse(zeros(257,randi([minInterpause maxInterpause],1,1)))};
end


% Now, generate the final sequence

for i = 1:numel(Sequenced)-1
    Sequenced2{i} = Sequenced{i};
end

% Now start compiling X
X = [];
for i = 1:size(Sequenced2,2)
    x = Sequenced2{i};
    X = [X x];
end

% Ok, now time to start compiling the full sng

Y = [];

for i = 1:size(X,2)
    Y = [Y X{i}];
end

Y = Y(:,1:180000); %clips to 3 minutes









