function [threshold,noiseMean,noisesd] = ...
    vocThreshold(wavfile,nfft,samprate,threshLength,freqRange,manual,pngPrint)
%% Read Samples
    s = wavread(wavfile,'size'); 
    samples = randsample(1+threshLength*samprate:s(1) - threshLength*samprate,1);
    samples = [samples, samples + threshLength*samprate - 1];    
    x = double(wavread(wavfile,samples,'native'));
    
%% Calculate the spectrogram of the block
    if size(x,2) < size(x,1)
        x = x';
    end
    padder = zeros(1,nfft/2);
    x = [padder,x];
    x = [x,padder(1:nfft/2-mod(length(x),nfft/2))];
    x1 = reshape(x,nfft/2,length(x)/(nfft/2));
    x1(:,end) = [];
    x2 = reshape(x,nfft/2,length(x)/(nfft/2));
    x2(:,1) = [];
    Y = [x1;x2];
    Z = fft(repmat(hamming(nfft),1,size(Y,2)).*Y);
    y = abs(Z((1:nfft/2+1),:));

    f = linspace(0,samprate/2,nfft/2 + 1);
    k{1} = find(f>freqRange(1),1);
    k{2} = find(f<freqRange(2));
    k{2} = k{2}(end);
    k = cat(2,k{:});
    
    y = y(k(1):k(2),:);
    y = log10(y);
    y = y(:);
    y = y(~isinf(y));
    
%% Automatic Threshold Setting

% Estimate mean and standard deviation of the noise, using the 68% rule and
% the histogram of the data. This assumes the peak is the noise. However it
% will be skewed because the distribution also contains real data. Thus we
% assume the noise is normally distributed and 68% of the data are
% contained within one standard deviation of the peak.

[h,b] = fdc(y','noplot');
h = h./sum(h);

[maxh,maxhi] = max(h);
hint = maxh;
c = 1;
while hint < 0.68
    hint = sum(h(maxhi-c:maxhi+c));
    c = c + 1;
end
noisesd = b(maxhi+c) - b(maxhi);
noiseMean = b(maxhi);
threshold = noiseMean + 3*noisesd;

if manual == 0
    if pngPrint == 1
        threshFigure = figure;
        plot(y); hold on;
        set(gca,'YLim',[-3 6]);
        plot([1,length(y)],[threshold,threshold],'r');
        print('-dpng',[wavfile,'_withthreshold.png']);
    end
%% Manual Threshold Setting
elseif manual == 1
    threshFigure = figure;
    plot(y);
    hold on;
    threshLine = plot([1,length(y)],[threshold,threshold],'r');
    
    ask = 2;
    
    while ask == 2
        w = warndlg('I`ve drawn a line. Inspect the graph and maybe pick a new value. If you like my line, just click OK and keep the default answer.');
        waitfor(w);
        newThresh = inputdlg({'Enter the threshold value'},'Spectrogram Threshold',1,{num2str(threshold)});
        delete(threshLine);
        
        threshold = str2double(newThresh);
        threshLine = plot([1,length(y)],[threshold,threshold],'r');
        
        askText = inputdlg({'I`ve re-drawn a line. Do you like it? (1) Yes or (2) No?'},'Spectrogram Threshold',1,{'1'});
        ask = str2double(askText);
    end   
    if pngPrint == 1
        set(gca,'YLim',[-3 6]);
        print('-dpng',[wavfile,'_withthreshold.png']);
    end
end
close(threshFigure);

threshold = 10.^threshold;
end
