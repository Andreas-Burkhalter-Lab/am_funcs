% fid = fopen('/mnt/fw1/d.bin','r');
% headersize = ReadAIHeader(fid);
% a = fread(fid, 'int16');
% plot (a);
% fclose(fid);


%  fid = fopen('/mnt/fw1/d.sng','r');
%  sng = ReadSonogram(fid);
%  absSng = abs(sng);
%  Bo = sngdecimate(absSng,2*2);
%  colormap(1-gray);
%  imagesc([0 10],[0 8192/2000],log10(Bo));
%  fclose(fid);

% fid = fopen('/mnt/fw1/v5.sng','r');
% sng = ReadSonogram(fid);
% absSng = abs(sng);
% figure
% plot(absSng);
% fclose(fid);

% fid = fopen('/mnt/fw1/sensordata.bin','r');
% a = fread(fid, 'float32');
% plot (a);
% fclose(fid);

% d = loadmc('~/m+f2.bin',[160 163]);
% plot(d);


%  fid = fopen('~/m+f2.bin','r');
%  headersize = ReadAIHeader(fid);
%  a = fread(fid, 'int16');
%  hfig = figure;
%  hax = axes('Parent',hfig);
%  x = linspace(0,180,length(a));
%  plot (x,a,'Parent',hax);
%  set(hax,'XLim',[160.85 161.2]);
%  fclose(fid);

% starts with low voltage
[transitions,initVoltage,header] = ReadSensor('/mnt/fw1/Y02-06-12/m6s13.det');
x = [0; transitions; 22282240];
len = ceil(length(x)/2);
Y = [zeros(1,len);ones(1,len)];
y = Y(1:length(x));
y(end)=y(end-1);
trans = (transitions/250000)'
secsperTransition = diff(transitions/250000)';

for (i=1:2:header.numTransitions-1)
   spt =  (secsperTransition(i))
    transitions(i)/250000   
    if (secsperTransition(i) > 0.02)
         c = secsperTransition(i)
         time = transitions(i)/250000
         break
    end
end
% index = find(secsperTransition > 1.45)
% transitions(index(1))/250000
[xb,yb] = stairs(x/250000,y);
plot(xb,yb);
sliderwindow(gca)

% starts with high voltage
% [transitions,initVoltage,header] = ReadSensor('m19sW.det');
% transitions(1)=[];
% x = [0; transitions; 22282240];
% len = ceil(length(x)/2);
% Y = [zeros(1,len);ones(1,len)];
% y = Y(1:length(x));
% y(end)=y(end-1);
% transitions/250000
% [xb,yb] = stairs(x/250000,y);
% plot(xb,yb);