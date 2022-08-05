function output = arriaga(pf,snips)

% AM 9/15/14: Note! This function analyzes mean frequency and standard
% deviation for all syllables types together, whereas Arriaga analyzed only
% mean frequency and standard deviation for type-A syllables. 
% I have edited this function to also calculate these variables for type-A
% syllables apart from other types of syllables. 

% Input pf is a 1xn cell array, where each cell is instantaneous peak frequency
% (bin with greatest FFT amplitude from sound2sng over times output by
% whistimes). This can be p.peakfreq (output of whisparams)

% Input snips is  1 x cell array with sparse sonograms for each syllable.
% Necessary for spectral purity calculations
% nfreq is number of frequency bins used in FFT, for spectral purity
% calculation.

% Note, Arriaga classes syllables based on freqrange of 35000 - 125000, and
% almost more importantly, based on a duration threshold of 10 ms. The
% whistimesdefaults uses a threshold of 5 ms, so only include syllables >=
% 10 ms duration if desired to fully reproduce Arriaga analysis.

% This function has been validated using real data from some of my
% recordings as well as some fake data with one of each type of syllable
% predesigned.

% The output from Tim's code should be truncated to remove any zeros or NaNs due to
% instantaneous loss of power in any time bin, but just in case:
count = 1;
for i = 1:length(pf)
    if ~isempty(pf{i})
    includes{count} = pf{i}(find(pf{i}~= 0 & ~isnan(pf{i})));
    count = count +1;
    end
end

% Arriaga categories A - K:
%   A: All jumps < 10 kHz
%   B: 1 negative jump, > 10 khz
%   C: 1 positive jump, > 10 khz
%   D: 2 negative jumps, > 10 khz
%   E: 1 negative followed by 1 positive, > 10 khz
%   F: 1 positive jump followed by 1 negative, > 10 khz
%   G: Negative, negative, positive
%   H: Positive, negative, positive, negative
%   I: ??? J: ??? K: ????     
% Classification vectors
% Algorithm:
%   class = jumps(jumps > 10000);
%   for i = 1: n elements of class
%       class(i) = class(i)/abs(class(i));
%   end
%   The resultant class should be of the form: [1,-1,...];
%   Then tests for inclusion into different classes match the following
%   categories (from Arriaga 2011 Thesis & PlosOne paper 2012 with jarvis)
%
%   Note: He classifies I,J,K syllables, but NOWHERE does he define them,
%   so I've not included them here. They will be grouped in other.

A = 0; 
% Type A matches SS from Holy analysis more or less. 
% Empirically, I actually think 10 Khz is a good cutoff for SS. 
% If you plot the distribution of abs(derivatives), you can see a sharp 
% fall-off at about 10 Khz and then a second peak. This is an alternative 
% to clustering the diagonal and jump clusters with manual ROIs.
B = -1;
C = 1;
D = [-1,-1];
E = [-1,1];
F = [1,-1];
G = [-1,-1,1];
H = [1,-1,1,-1];

types = [1,2,3,4,5,6,7,8,9];
typenames = {'a','b','c','d','e','f','g','h','o'};

% Arriaga analysis:
%   1) Classification into A-K (but really I can only do A-H).
%   2) Calculation of fraction of each syllable type. We'll output both n
%   of each type & n./sum(n) for fractions.
%   3) Time averaged peak frequency (mean frequency) per syllable
%   4) Standard deviation of peak frequencies per syllable
%   5) For plotting: distribution of instantaneous peak frequencies per
%   syllable as well as distribution over each syllable class (lumped all
%   together it seems).
%   6) Starting frequency: frequency at bins(1)
%   7) 1/4 frequency: Frequency at bins(1+floor(length(bins)/4)) per syllable
%   8) 3/4 frequency: Frequency at bins(1+3*floor(length(bins)/4)) per syllable
%   9) Final frequency: frequency at bins(end)
%   10) Min instantaneous frequency
%   11) Max instantaneous frequency
%   12) Frequency variance = std(bins)^2
%   13) Spectral purity: fraction of all frequencies (from sngparms.nfreq)
%    used. Can calculate as length(unique(bins))/nfreq. This is the total
%    purity. I'm not sure however, whether he means the time averaged 
%    purity. Could calculate both.
%   

% Initialize storage for output
output.class = zeros(1,length(includes)); % from types
output.names = cell(1,length(includes));  % text label from typenames
output.classN = zeros(1,9); %Calculates number for 1-8 & other
output.classFraction = zeros(1,9);
output.classfreqDist = cell(1,9); % FDC distibution of frequencies by type
output.syllfreqDist = cell(1,length(includes)); % FDC distribution of frequencies by syllable
output.freqDist = []; % FDC total distribution over all syllables
output.means = zeros(1,length(includes)); % Time averaged frequency
output.mins = zeros(1,length(includes)); % min value for fundamental
output.maxs = zeros(1,length(includes)); % max value for fundamental
output.bandwidth = zeros(1,length(includes)); % range of fundamental
output.sds_all = zeros(1,length(includes)); % Standard deviation
output.var_all = zeros(1,length(includes)); % Variance
output.start = zeros(1,length(includes)); % t(1)
output.quarter1 = zeros(1,length(includes)); %t(1+floor(nbins/4)
output.quarter3 = zeros(1,length(includes)); %t(1+3*floor(nbins/4)
output.final = zeros(1,length(includes)); % t(end)
output.sp = zeros(1,length(includes)); % Total spectral purity by whistles
output.spav = zeros(1,length(includes)); % Average spectral purity by whistle


%Ok! we're ready to get started. For QC:
% Stored variable under arriaga/arriagatest.mat has two arrays:
% arriagatest:  Fake data
%   Syllable 1: Type A
%   Syllable 2: Type A, has some zeros
%   Syllable 3: Type A, has some NaNs
%   Syllable 4: Type B
%   Syllable 5: Type C
%   Syllable 6: Type D
%   Syllable 7: Type E
%   Syllable 8: Type F
%   Syllable 9: Type G
%   Syllable 10: Type H
%   Syllable 11: Other
% arriagatest2: Real data from a WT Bl6 Male/Female interaction from the
% recent Celf6 data set. Data was a sexually naive WT C57BL6/J adult male 
% paired with WT C57Bl6/J female, a few bouts from about 30 seconds of
% vocalization.
% Other variables are the whistle snips, twhis, and whisparams outputs for the
% test.sng file, and a version downsampled to 10000/sec.
% 



% Classify the syllables

%% Calculate the classifier vectors for each whistle:
for i = 1:length(includes)
    classifier = diff(includes{i}); % calculate derivative
    classifier = classifier(abs(classifier) >= 10000); % find >10kHz jumps
    classifier = classifier./abs(classifier); % normalize
    
    % classifier should be of the form A-H above. Test for equality:
    if isempty(classifier)
        classifier = 0;
    end
    % Apparently empty is somehow different from [], so you can't perform equality
    % tests.
    if numel(classifier) == 1 && classifier == A
        output.class(i) = types(1);
        output.names{i} = typenames{1};
    elseif numel(classifier) == 1 && classifier == B
        output.class(i) = types(2);
        output.names{i} = typenames{2};
    elseif numel(classifier) == 1 && classifier == C
        output.class(i) = types(3);
        output.names{i} = typenames{3};
    elseif numel(classifier) == 2 && sum(classifier == D) == 2
        output.class(i) = types(4);
        output.names{i} = typenames{4};
    elseif numel(classifier) == 2 && sum(classifier == E) == 2       
        output.class(i) = types(5);
        output.names{i} = typenames{5};
    elseif numel(classifier) == 2 && sum(classifier == F) == 2 
        output.class(i) = types(6);
        output.names{i} = typenames{6};
    elseif numel(classifier) == 3 && sum(classifier == G) == 3 
        output.class(i) = types(7);
        output.names{i} = typenames{7};
    elseif numel(classifier) == 4 && sum(classifier == H) == 4
        output.class(i) = types(8);
        output.names{i} = typenames{8};
    else
        output.class(i) = types(9);
        output.names{i} = typenames{9};
    end
 
    clear classifier;
end

%% Create cell array containing only type A syllables (AM added 9/15/14). 
includes_typeA = includes(find(cell2mat(output.names) == 'a'));

%% Calculate class fractions
for i = 1:9
output.classN(i) = numel(find(output.class == i));
end
output.classFraction = (output.classN)./sum(output.classN);

%% Calculate frequency distributions by syllable type
bins = linspace(0,125000,85); % This might not really be the right number of bins but it's good for now and standardizes all the distributions.
% It's based roughly on friedman diaconis bin size I've used with real
% data. binsize ~1500 Hz

for i = 1:length(types)
    typeData = includes(output.class == types(i));
    typeData = cat(2,typeData{:});
    output.classfreqDist{i} = hist(typeData,bins);
    output.classfreqDist{i} = (output.classfreqDist{i})./sum(output.classfreqDist{i}); 
    clear typeData;
end

%% Calculate frequency distributions by syllable
for i = 1:length(includes)
    output.syllfreqDist{i} = hist(includes{i},bins);
    output.syllfreqDist{i} = (output.syllfreqDist{i})./sum(output.syllfreqDist{i}); 
    clear typeData;
end


%% Calculate total frequency distribution

output.freqDist = hist(cat(2,includes{:}),bins);
output.freqDist = (output.freqDist)./sum(output.freqDist);

%% Calculate rest of output:

% Note! For spectral purity, I didn't make a fake sparse sonogram to match
% arriagatest (the variable with fake peak freqs). For illustration, use
% arriagasnips2 (matches arriagatest2 data)

for i = 1:length(includes)
    output.means(i) = mean(includes{i});
    output.mins(i) = mean(min(includes{i}));
    output.maxs(i) = mean(max(includes{i}));
    output.bandwidth(i) = output.maxs(i) - output.mins(i);
    output.sds_all(i) = std(includes{i});                       %% AM changed 9/15/14
    output.var_all(i) = output.sds_all(i)^2;                       %% AM changed 9/15/14
    output.start(i) = includes{i}(1);
    output.quarter1(i) = includes{i}(1 + floor(length(includes{i})/4));
    output.quarter3(i) = includes{i}(1 + 3*floor(length(includes{i})/4));
    output.final(i) = includes{i}(end);
    % For SP calculations:
    powerspec = abs(snips{i}).^2;
    % Truncate any instantaneous power loss and average, divided by nfreq:\
    totpow = sum(full(powerspec));
    maxpow = max(full(powerspec));
    output.spav(i) = mean(maxpow(find(totpow))./totpow(find(totpow)));
    % Find the total specpurity over the whole syllable. This isn't just sum(numFreqsAccum) because the same freqs will be used from column to column;
    % First integrate the power, then calculate:
    output.sp(i) = max(sum(full(powerspec),2))/sum(sum(full(powerspec),2));
    clear powerspec totpow maxpow;
    % spav is better - it's the average spectral purity over time. It's how
    % Tim's code calculates whether a time bin belongs in the syllable or
    % not. sp tells you how many frequency bins overall the syllable
    % employs. But again, I don't know what Arriaga means so I've
    % calculated both. sp may be >0.25, but spav is never greater than 0.25
    % because that's the threshold in whistimesdefaults.
end

%% Loop added by AM 9/15/14 to calculate mean frequency and intra-syllable
%  standard deviation of type-A syllables apart from other syllable types. 
for i = 1:length(includes_typeA)
    output.sds_typeA(i) = std(includes_typeA{i});
    output.var_typeA(i) = output.sds_typeA(i)^2;  
end

%% Here is an output of classification using arriagatest2 real data:

% His mean frequency for type A in C57 animals: ~75 kHz. Mine from the
% sample of 34 type A syllables (SS) 76.5 kHz

% % classFraction: [0.6182 0.1273 0.0727 0 0.0727 0.0182 0 0 0.0909]

% The expected values from his thesis are given as a pie chart without
% ennumeration. But at least for A,B,E from BxD sham surgery animals:

% Type A = ~65%
% Type B = ~10%
% Type E = ~10%

% My Type A & B are very close from the 55 syllable test data. E is a
% little lower but probably in the range of whatever variability he sees.
%
% I have a about 10% in this small sample of unclassified by his method.
%

% The means of his type A syllable standard deviations is about 10 kHz. 
% From my sample of 55 syllables, 34 were type A, and most of mine had a
% standard deviation of only about 2 kHz.
%

% Overall, I think this probably reproduces his analysis well, so I can go
% ahead and run it on a larger data set.



