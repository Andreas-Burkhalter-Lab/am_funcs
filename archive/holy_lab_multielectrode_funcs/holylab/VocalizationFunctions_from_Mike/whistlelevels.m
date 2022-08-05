function [avamp,dbavamp,ref] = whistlelevels(twhisfilename,twhisdir,wavfilename,wavdir,nfreq,samprate)
%% These are inputs to the function - to be removed later.
df = (samprate/2)/nfreq;
sngparms = struct('plot',0,'nfreq',nfreq,'freqrange',[0 (samprate/2)],...
                    'progressbar',0,'threshold',0);
tempdirname = 'Temp';
    if ~isdir(tempdirname);
        mkdir(tempdirname);
    end
%% Load twhis

load(['./',twhisdir,'/',twhisfilename]); 

% Stores start & end times calculated for whistles as variable 'twhis'

if size(twhis,2) >= 1

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
fclose('all'); % close any open files.
whistlesngs = dir([strrep(wavfilename,'.WAV',''),'*.SNG']);

%% Calculate the time-averaged amplitude per frequency bin for each sonogram

avamp = cell(1,numel(whistlesngs));
for i = 1:numel(whistlesngs)
    [sng,header] = ReadSonogram(whistlesngs(i).name);
    df = ((header.scanrate/2)/(header.nfreq-1));
    sng = abs(sng)';
    avamp{i} = mean(full(sng));
end

avamp = cat(1,avamp{:});
mins = min(avamp);
if length(mins) > 1
ref = min(mins(floor(25000/df):floor(110000/df)));
else
ref = mins;
end

% For dB calculation:
%   reference amplitude is the minimum time averaged amplitude from among
%   all the whistles in the whistle portion of the spectrum based on Tim's paper
%   25000 - 110000 Hz.

dbavamp = 20*log10(avamp./ref);

cd ..;

else
    avamp = zeros(1,nfreq+1);
    avamp(1:end) = NaN;
    dbavamp = zeros(1,nfreq+1);
    dbavamp(1:end) = NaN;
    ref = NaN;
end

end







