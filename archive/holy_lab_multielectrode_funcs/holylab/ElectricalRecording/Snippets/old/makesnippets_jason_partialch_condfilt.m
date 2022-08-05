merecfiles = dirbytime('*.merec_jason');
merecfiles = merecfiles(1:end);
[header,fid] = readheader(merecfiles{1});
w = fread(fid,[header.numch 50000],'*int16');

so = struct('polarity', -1);
so.Fs = header.scanrate; % Fs is scanrate
so.Hz60 = 1; % need set this, what in old snippetfileauto.m is incomplete
so = snipoptions(so); % this is not necessary because in snippetfile.m it is also called
chanIndx = 1:size(header.channels,2);
w = filterint16(so.condfiltb,so.condfilta,w,chanIndx);
w=double(w);

mn = median(w,2);
wf = w - repmat(mn,1,size(w,2));
mna = mean(abs(wf),2);
threshL = -6*mna + mn;
threshU =  6*mna + mn;
% thresh = thresh';
thresh = [threshL'; threshU'];
%tChannels=header.channels;
tChannels=[25 33];
[tComCh, tIndex]=intersect(header.channels, tChannels);
tIndex=sort(tIndex);
thresh=thresh(:, tIndex);

tScanRange=[3000 header.nscans-1];

for i = 1:length(merecfiles)
  [pathstr,basename,ext] = fileparts(merecfiles{i});
  basename
  snippetfile(merecfiles{i}, tScanRange,tChannels,thresh,[-10 30],[basename,'.ssnp'],so);
end
