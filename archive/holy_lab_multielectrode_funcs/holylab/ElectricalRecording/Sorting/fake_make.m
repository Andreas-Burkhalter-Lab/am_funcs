sniprange = [-10 30];
sniplen = diff(sniprange)+1;
x = linspace(0,2*pi,sniplen)';
snipclass = [sin(x).*exp(-x),-sin(x).*exp(-x),2*sin(x).*exp(-x),...
             sin(x).^2.*exp(-x),zeros(size(x))];
nclasses = size(snipclass,2);
figure
for i = 1:nclasses
  col = unique_color(i,nclasses);
  line(sniprange(1):sniprange(end),snipclass(:,i),'Color',col)
end

rate = rand(1,nclasses)+0.2;
minspikes = 200;
duration = minspikes/min(rate);
nspikes = round(duration*rate);
noiseamp = 0.025;
nfiles = 2;
noutliers = 3;
outlieramp = 0.2;
channel = 55;

load fake_baseheader  % This loads the variable "h"
tstart = now;

for i = 1:nfiles
  figure
  t = {};
  snip = {};
  clust = {};
  for j = 1:nclasses
    t{j} = round(h.scanrate*duration*rand(1,nspikes(j)));
    snip{j} = repmat(snipclass(:,j),1,nspikes(j)) + ...
              noiseamp*randn(sniplen,nspikes(j));
    snip{j} = round((snip{j} - h.scaleoff)/h.scalemult);
    clust{j} = repmat(j,1,nspikes(j));
    col = unique_color(j,nclasses);
    line(sniprange(1):sniprange(end),snip{j},'Color',col)
  end
  % Add some outliers
  t{end+1} = round(h.scanrate*duration*rand(1,noutliers));
  snip{end+1} = round((outlieramp*randn(sniplen,noutliers)-h.scaleoff)/h.scalemult);
  clust{end+1} = nan(1,noutliers);
  % Now sort in temporal order
  tall = cat(2,t{:});
  [tall,swapIndex] = sort(tall);
  snipall = cat(2,snip{:});
  clustall = cat(2,clust{:});
  snipall = snipall(:,swapIndex);
  clustall = clustall(:,swapIndex);
  clust_in_file{i} = clustall;
  finetimes = zeros(size(tall));
  thresh = -h.scaleoff/h.scalemult;
  [detpeaks,dpIndex] = max(abs(snipall - thresh));
  idx = sub2ind(size(snipall),dpIndex,1:size(snipall,2));
  detpeaks = snipall(idx);
  % Get the header set up
  tHead = h.wholeheader;
  tHead = update_value(tHead,'nscans', num2str(ceil(duration*h.scanrate)));
  tHead = update_value(tHead,'datetime', ...
                       datestr(tstart + 1.1*duration*(i-1)/(24*3600)));
  tHead = update_value(tHead,'snipbeginoffset', num2str(sniprange(1)));
  tHead = update_value(tHead,'snipendoffset', num2str(sniprange(2)));
  tHead = update_value(tHead,'thresh', num2str(thresh([1 1])));
  tHead = update_value(tHead,'polarity', num2str(0));
  tHead = update_value(tHead,'snippet input file', ...
                       ['fake' num2str(i) '.merec']);
  % Write the snippet file
  fname = ['fake' num2str(i) '.ssnp'];
  fid = fopen(fname,'w');
  if (fid < 0)
    error(['Error opening fake snippet file ' fname]);
  end
  update_header(fid,tHead);
  snipfile_append_channel(fid,tHead, channel, tall, snipall, finetimes, ...
                          detpeaks);
  fclose(fid);
end    
current_sort_info = struct('channel',channel,'sniprange',sniprange,...
                           'use_projection',0,'t2V',0,...
                           'landmarkWaveform',snipclass,...
                           'landmarkT',zeros(1,nclasses));
save 'clust_rightanswer' clust_in_file current_sort_info

