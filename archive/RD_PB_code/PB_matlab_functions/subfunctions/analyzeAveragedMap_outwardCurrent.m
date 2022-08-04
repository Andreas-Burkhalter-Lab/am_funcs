function map=analyzeAveragedMap_outwardCurrent(map)

% last edited by AM 2/19/16 on recording comp
% script now saves full grid and mean of peak currents and charge integral

% Analyses a single map (usually a mean maP and returns arrays with the
% results and plots of traces and the analysis.
%   LTP 2006
% MH 20100406

[rows, numTraces] = size(map.averageMap);
sr=map.samplingRate;

%2009-11-05 for filtering AP
filter_flag = 0;freq = 100;sampleRate = 10000;
if filter_flag
    [b,a] = butter(4,freq*2/sampleRate,'low');
    Hd = dfilt.df2t(b,a);
    for i = 1:size(map.averageMap,2)
        map.averageMap(:,i) = filter(Hd,map.averageMap(:,i))';
    end
end

baselineStartIndex = 1;
baselineEndIndex = 999;

%SET: STIMULUS 'ON' TIME
%SET: OTHER PARAMETERS FOR RESPONSE ONSET AND SYN RESPONSE

%stimOnInd=1000; %<--------------1000!
%responseOnsetEnd=1230;
%synEndInd=1750;
%chargeEndInd=2500;

stimOnInd=1000; %<--------------1000!
responseOnsetEnd=1150;
synEndInd=1750;
chargeEndInd=2500;

baselineMedians = median(map.averageMap(1:999,:));
baselineSDs = std(map.averageMap(baselineStartIndex:baselineEndIndex,:));
SD=mean(baselineSDs);

% THRESHOLDS ****************************************************************
% NB: these are vectors

dirLevel=3; %<---------Config

dirNegThresh = baselineMedians+dirLevel*SD; %% AM changed -dirlevel to +dirLevel
threshold=repmat(dirNegThresh, rows, 1);
for i=1:numTraces
    % AM changed < to > in following if statement
    if isempty(find(diff(map.averageMap(stimOnInd:responseOnsetEnd,i)>threshold(stimOnInd:responseOnsetEnd,i))==1));
        onset(i)=NaN;
        Maximum(i)=max(map.averageMap(stimOnInd:synEndInd,i));
        minimum(i)=NaN;
        minOnset(i)=NaN;
        onset90(i)=NaN;
        onset10(i)=NaN;
        riseTime1090(i)=NaN;
        traceMean(i)=0; % for printing the significant traces TM, comment out if need the whole map
    else
        % AM changed < to > in following line; min finds first instance of
        % slope threshold crossing
        onset(i)=min(find(diff(map.averageMap(stimOnInd:responseOnsetEnd,i)>threshold(stimOnInd:responseOnsetEnd,i))==1))/sr;
        [Minimum MinOnset]=min(map.averageMap(stimOnInd:synEndInd,i));
        Maximum(i)=max(map.averageMap(stimOnInd:synEndInd,i));

        % if singleEventEPSCdiscriminator(map.averageMap(stimOnInd:synEndInd,i),onset(i)*sr)  %only analyse riseTime for single component event;
        %     onset90(i)=NaN;
        %     onset10(i)=NaN;
        % else
        
        % AM changed < to > and min((map. to max((map. in following line; this line finds the first
        % instance of trace exceeding 90% of max voltage
        onset90(i)= min(find(diff((map.averageMap(stimOnInd:synEndInd,i)>.9*max((map.averageMap(stimOnInd:synEndInd,i))))==1)));
        try
            % AM changed < to > and min((map. to max((map. in following line; this line finds the first
            % instance of trace exceeding 10% of max voltage
            onset10(i)=min(find(diff((map.averageMap(stimOnInd:synEndInd,i)>.1*max((map.averageMap(stimOnInd:synEndInd,i))))==1)));
        catch
            onset10(i)=NaN;
        end
        % end
        minimum(i)=Minimum;
        minOnset(i)=MinOnset/sr;
        riseTime1090(i)=(onset90(i)-onset10(i))/sr;
        traceMean(i)=mean(map.averageMap(stimOnInd:synEndInd,i)); % comment out it need the whole map TM
    end
%     traceMean(i)=mean(map.averageMap(stimOnInd:synEndInd,i));
    integral(i)= trapz(map.averageMap(stimOnInd:chargeEndInd,i))/sr;

    [a b]=find(map.pattern{1}==i);
    if a==11&&b==8
        iiii =  0;
    end

    mapOnset(a,b)=onset(i);
    mapMinOnset(a,b)=minOnset(i);
    mapMin(a,b)=minimum(i);
    mapMean(a,b)=traceMean(i);
    mapIntegral(a,b)=integral(i);
    mapRiseTime1090(a,b)=riseTime1090(i);
    mapMaximum(a,b)=Maximum(i);
end
map.analysis.xSpacing=map.xSpacing;
map.analysis.ySpacing=map.ySpacing;
map.analysis.onset=mapOnset;
map.analysis.minOnset=mapMinOnset;
map.analysis.minimum=mapMin;
map.analysis.mean=mapMean;
map.analysis.integral=mapIntegral;
map.analysis.riseTime1090=mapRiseTime1090;
map.analysis.maximum=mapMaximum;
%[map.soma1Coordinates(1)
%map.soma1Coordinates(2)]=somaPositionTransformerMH(map.soma1Coordinates(1)
%, map.soma1Coordinates(2), map.spatialRotation,0,0) edited by MH to plot
%correctly
[map.soma1Coordinates(1) map.soma1Coordinates(2)]=somaPositionTransformerMH(map.soma1Coordinates(1), map.soma1Coordinates(2), map.spatialRotation,map.xPatternOffset,map.yPatternOffset);

%Address issue of different size maps: some have 12 columns, others have 8
%or 10.  When <12, add columns of NaN to make the same shape MH20100412

for n=1:size(map.analysis.mean,1)
    nancolumnvector(n)=NaN;
end

% if size(map.analysis.mean,2)==8
%     map.analysis.mean=horzcat(nancolumnvector',nancolumnvector',map.analysis.mean,nancolumnvector',nancolumnvector');
%     map.analysis.minimum=horzcat(nancolumnvector',nancolumnvector',map.analysis.minimum,nancolumnvector',nancolumnvector');
%     map.analysis.integral=horzcat(nancolumnvector',nancolumnvector',map.analysis.integral,nancolumnvector',nancolumnvector');
%     map.analysis.onset=horzcat(nancolumnvector',nancolumnvector',map.analysis.onset,nancolumnvector',nancolumnvector');
%     map.analysis.minOnset=horzcat(nancolumnvector',nancolumnvector',map.analysis.minOnset,nancolumnvector',nancolumnvector');
%     map.analysis.riseTime1090=horzcat(nancolumnvector',nancolumnvector',map.analysis.riseTime1090,nancolumnvector',nancolumnvector');
% end
% 
% if size(map.analysis.mean,2)==10
%     map.analysis.mean=horzcat(nancolumnvector',map.analysis.mean,nancolumnvector');
%     map.analysis.minimum=horzcat(nancolumnvector',map.analysis.minimum,nancolumnvector');
%     map.analysis.integral=horzcat(nancolumnvector',map.analysis.integral,nancolumnvector');
%     map.analysis.onset=horzcat(nancolumnvector',map.analysis.onset,nancolumnvector');
%     map.analysis.minOnset=horzcat(nancolumnvector',map.analysis.minOnset,nancolumnvector');
%     map.analysis.riseTime1090=horzcat(nancolumnvector',map.analysis.riseTime1090,nancolumnvector');
% end


%--------------display Analysis----------------------
% clc
%
% disp('function = TM');
% disp(['map.experimentNumber= ',map.experimentNumber]);
% disp('map.presynaptic=');
% disp('map.postsynaptic=');
% disp('%--------------------------------------------------------');
% disp('map.distanceToPia=');
% disp('% map.distanceToBarrelTop/Bottom=');
% disp('map.cortexThickness=');
% disp('map.layer1Row=');
% disp('% map.distanceToMidline=');
% disp('% map.GCCToMidline=');
% disp('%--------------------------------------------------------');
% fprintf('map.soma1Coordinates=[ %g %g ]; \n', [map.soma1Coordinates(1)-map.xPatternOffset map.soma1Coordinates(2)-map.yPatternOffset])
% fprintf('map.laserPower= %g ; \n', mean(map.laserPower))
% disp('map.pulseWidth=');
% fprintf('map.xSpacing= %g ; \n', map.xSpacing)
% fprintf('map.ySpacing= %g ; \n', map.ySpacing)
% fprintf('map.numberOfMaps= %g ; \n', length(map.laserPower));
% fprintf('map.analysisWindow= %g ; \n', (synEndInd-stimOnInd)/sr);
%
%
% disp('map.mean=[')
% disp(map.analysis.mean)
% disp('];')
%
% disp('map.min=[')
% disp(map.analysis.minimum)
% disp('];')
%
% disp('map.integral=[')
% disp(map.analysis.integral)
% disp('];')
%
% disp('map.onset=[')
% disp(map.analysis.onset)
% disp('];')
%
% disp('map.minOnset=[')
% disp(map.analysis.minOnset)
% disp('];')
%
% disp('map.riseTime1090=[')
% disp(map.analysis.riseTime1090)
% disp('];')

%%% By HN for directly input into file

%[filename, pathname] = uiputfile(['C:\hooksm\_Cartography\dataMfiles\',map.experimentNumber,'.m'], 'Save data file as');
[filename, pathname] = uiputfile(['E:\HHMI 2011 Oct23 - Dec22\Analysis\Weiguo Data\DataMfiles\',map.experimentNumber,'.m'], 'Save data file as');

fid = fopen(fullfile(pathname, filename),'wt');
fprintf(fid, 'function map = %s\n', map.experimentNumber);
fprintf(fid, '%s\n','%Maps\n');
fprintf(fid, 'map.experimentNumber= \''%s\'';\n', map.experimentNumber);
fprintf(fid, 'map.presynaptic='''';\n');
fprintf(fid, 'map.postsynaptic='''';\n');
fprintf(fid, 'map.hemisphere='''';\n');
fprintf(fid, '%s\n','%--------------------------------------------------------');
fprintf(fid, 'map.distanceToPia= nan;\n');
fprintf(fid, '%s\n','% map.distanceToBarrelTop/Bottom=');
fprintf(fid, 'map.cortexThickness= nan;\n');
fprintf(fid, 'map.yfrac= nan;\n');
fprintf(fid, 'map.layer1Row= nan;\n');
fprintf(fid, '%s\n','% map.distanceToMidline=');
fprintf(fid, '%s\n','% map.GCCToMidline=');
fprintf(fid, '%s\n','%--------------------------------------------------------');
%fprintf(fid, 'map.soma1Coordinates=[ %g %g ]; \n', [map.soma1Coordinates(1)-map.xPatternOffset...
%    map.soma1Coordinates(2)-map.yPatternOffset]); edited to correct soma
%    coordinates for plotting
fprintf(fid, 'map.soma1Coordinates=[ %g %g ]; \n', [map.soma1Coordinates(1)...
    map.soma1Coordinates(2)]);
fprintf(fid, 'map.laserPower= %g ; \n', mean(map.laserPower));
fprintf(fid, 'map.pulseWidth= nan;\n');
fprintf(fid, 'map.xSpacing= %g ; \n', map.xSpacing);
fprintf(fid, 'map.ySpacing= %g ; \n', map.ySpacing);
fprintf(fid, 'map.numberOfMaps= %g ; \n', length(map.laserPower));
fprintf(fid, 'map.analysisWindow= %g ; \n', (synEndInd-stimOnInd)/sr);

fprintf(fid, ['map.mean=[']);
for n=1:size(map.analysis.mean,1)
    fprintf(fid, [num2str(map.analysis.mean(n,:))]);
    fprintf(fid, [';\n']);
end
fprintf(fid, ['];']);
fprintf(fid, ['\n']);

fprintf(fid, ['map.min=[']);
for n=1:size(map.analysis.minimum,1)
    fprintf(fid, [num2str(map.analysis.minimum(n,:))]);
    fprintf(fid, [';\n']);
end
fprintf(fid, ['];']);
fprintf(fid, ['\n']);

fprintf(fid, ['map.max=[']);
for n=1:size(map.analysis.maximum,1)
    fprintf(fid, [num2str(map.analysis.maximum(n,:))]);
    fprintf(fid, [';\n']);
end
fprintf(fid, ['];']);
fprintf(fid, ['\n']);

fprintf(fid, ['map.integral=[']);
for n=1:size(map.analysis.integral,1)
    fprintf(fid, [num2str(map.analysis.integral(n,:))]);
    fprintf(fid, [';\n']);
end
fprintf(fid, ['];']);
fprintf(fid, ['\n']);

fprintf(fid, ['map.onset=[']);
for n=1:size(map.analysis.onset,1)
    fprintf(fid, [num2str(map.analysis.onset(n,:))]);
    fprintf(fid, [';\n']);
end
fprintf(fid, ['];']);
fprintf(fid, ['\n']);

fprintf(fid, ['map.minOnset=[']);
for n=1:size(map.analysis.minOnset,1)
    fprintf(fid, [num2str(map.analysis.minOnset(n,:))]);
    fprintf(fid, [';\n']);
end
fprintf(fid, ['];']);
fprintf(fid, ['\n']);

fprintf(fid, ['map.riseTime1090=[']);
for n=1:size(map.analysis.riseTime1090,1)
    fprintf(fid, [num2str(map.analysis.riseTime1090(n,:))]);
    fprintf(fid, [';\n']);
end
fprintf(fid, ['];']);
fprintf(fid, ['\n']);

fclose(fid);

%%% MH changed 20100421
%%% By HN for directly input into file

mapAnalysisPlotTTXMH(map)
arrayTracesAsInputMapTTX(map)

% section added by AM
peaks_grid = map.analysis.maximum;
charge_integral_grid = map.analysis.integral;
peak_mean = mean(mean(peaks_grid));
charge_integral_mean = mean(mean(charge_integral_grid));
save('analysis_results','peaks_grid','charge_integral_grid','peak_mean','charge_integral_mean');


