%% These are inputs to the function - to be removed later.

twhisfilename = 'twhisPET101_190.WAV.SNG.mat';
twhisdir = 'Outputs';
wavfilename = 'PET101_190.WAV';
wavdir = 'WAVs';
sngparms = struct('plot',0,'nfreq',nfreq,'freqrange',[0 (samprate/2)],...
                    'progressbar',0,'threshold',0);
tempdirname = 'Temp';
  
    if ~isdir(tempdirname);
        mkdir(tempdirname);
    end
    
    
samprate = 250000;
nfreq = 256;
df = (samprate/2)/nfreq;


%% Load twhis

load(['./',twhisdir,'/',twhisfilename]); % Stores start & end times calculated for whistles as variable 'twhis'
whissamples = zeros(1,size(twhis,2));
whissamples(1,:) = floor(twhis(1,:)*samprate);
whissamples(2,:) = ceil(twhis(2,:)*samprate);

%% Make unthresholded sonograms

for i = 1:size(whissamples,2)
    [y,fs,nbits] = wavread(['./',wavdir,'/',wavfilename],[whissamples(1,i),whissamples(2,i)]);
    wavwrite(y,fs,nbits,['./',tempdirname,'/',strrep(wavfilename,'.WAV',''),'.',sprintf('%05d',i),'.WAV']);
end

cd(tempdirname);
whistleclips = dir([strrep(wavfilename,'.WAV',''),'*.WAV']);
for i = 1:numel(whistleclips)
    sound2sng(whistleclips(i).name,sngparms,[strrep(whistleclips(i).name,'.WAV',''),'.SNG']);
end

whistlesngs = dir([strrep(wavfilename,'.WAV',''),'*.SNG']);

%% 
% Calculate the :
%   avamp: time-averaged amplitude per frequency bin for each sonogram
%   totamp: time-integrated amplitude per frequency bin for each sonogram
%   mins: reference values for min amplitude
%   

avamp = cell(1,numel(whistlesngs));
for i = 1:numel(whistlesngs)
    [sng,header] = ReadSonogram(whistlesngs(i).name);
    df = ((header.scanrate/2)/(header.nfreq-1));
    sng = abs(sng)';
    avamp{i} = mean(full(sng));
end

avamp = cat(1,avamp{:});
mins = min(avamp);
ref = min(mins(floor(25000/df):floor(110000/df)));
% For dB calculation:
%   reference amplitude is the minimum time averaged amplitude from among
%   all the whistles in the whistle portion of the spectrum 25000 - 110000

dbavamp = 20*log10(avamp./ref);







