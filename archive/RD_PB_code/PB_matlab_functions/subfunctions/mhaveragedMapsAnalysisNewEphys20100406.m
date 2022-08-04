function map=mhaveragedMapsAnalysisNewEphys20100406

%Loads a number of maps from the same cell and averages the maps.
%Call this function to start the analysis of an experiment.
% Leopoldo Petreanu 2006
% KS011708 fixed a bug with averaging; output the traceStack
% MH 20100406

numberOfMaps=input('number of maps?');
for j=1:numberOfMaps;
    loadedData{j}=loadMapNewEphysMH; %load data
    laserPower(j)=loadedData{j}.laserPower;
    [a,b]=size(loadedData{j}.pattern); %uncaging pattern
    patternRotation=str2num(loadedData{j}.patternRotation);
    patternFlip=loadedData{j}.patternFlip;
    pattern{j}=loadedData{j}.pattern;
    % % if patternFlip
    % %     pattern{j} = flipdim(pattern{j}, 2);
    % % else
    % switch patternRotation
    %     case '0'
    %         %Do nothing.
    %     case '90'
    %         pattern{j} = rot90(pattern{j}, 1);
    %     case '180'
    %         pattern{j} = rot90(pattern{j}, 2);
    %     case '270'
    %         pattern = rot90(pattern{j}, 3);
    % end
    % %end
    meanData=zeros(a,b);
    [r c]=size(loadedData{j}.traces);
    dataArray{j}=loadedData{j}.traces;
    [rows, numTraces] = size(dataArray{j});
end

% BASELINE settings, mean, sd, baseline-subtracted array
baselineStartIndex = 1;
baselineEndIndex = 999;

for j=1:numberOfMaps
    baselineMeans = median(dataArray{j}(baselineStartIndex:baselineEndIndex,:));
    baselineSDs = std(dataArray{j}(baselineStartIndex:baselineEndIndex,:));
    SD=mean(baselineSDs);
    % baseline-subtracted array:
    baselineArray = repmat(baselineMeans, rows, 1);
    bsArray{j} = dataArray{j} - baselineArray;
end

% switch numberOfMaps
%     case 3
%         for i=1:numTraces
%    [a b]=find(pattern{1}==i);
%    traceMap2=pattern{2}(a,b);
%    traceMap3=pattern{3}(a,b);
%    traceStack(:,1)=bsArray{1}(:,i);
%    traceStack(:,2)=bsArray{2}(:,traceMap2);
%    traceStack(:,3)=bsArray{3}(:,traceMap3);
%    averageMap(:,i)=mean(traceStack,2);
%         end
%
%     case 2
%         for i=1:numTraces
%       [a b]=find(pattern{1}==i);
%    traceMap2=pattern{2}(a,b);
%
%    traceStack(:,1)=bsArray{1}(:,i);
%    traceStack(:,2)=bsArray{2}(:,traceMap2);
%
%    averageMap(:,i)=mean(traceStack,2);
%  end

if numberOfMaps == 1
    traceStack=zeros([size(bsArray{1}), numberOfMaps]); % KS011708 initialize the traceStack
    traceStack(:,:,1)=bsArray{1};
    averageMap = bsArray{1};
else
    traceStack=zeros([size(bsArray{1}), numberOfMaps]); % KS011708 initialize the traceStack
    traceStack(:,:,1)=bsArray{1};
    for j=2:numberOfMaps
        for i=1:numTraces
            [a b]=find(pattern{1}==i);
            trace=bsArray{j}(:,pattern{j}(a,b));
            traceStack(:,i,j)=trace;
        end
%        whos traceStack
        averageMap=mean(traceStack,3); 
%        stdMap=std(traceStack,3);
    end
end

%pattern1=pattern{1};
map.pattern=pattern;
map.soma1Coordinates=loadedData{1}.soma1Coordinates;
map.xSpacing=loadedData{1}.xSpacing;
map.ySpacing=loadedData{1}.ySpacing;
map.xPatternOffset=loadedData{1}.xPatternOffset;
map.yPatternOffset=loadedData{1}.yPatternOffset;
map.laserPower=laserPower;
map.traces=dataArray;
map.averageMap=averageMap;
map.maps=traceStack; % KS011708 save individual maps
map.samplingRate=loadedData{1}.samplingRate;
map.experimentNumber=loadedData{1}.experimentNumber;
map.spatialRotation=loadedData{1}.spatialRotation;
map.horizontalVideoDistance=loadedData{1}.horizontalVideoDistance;
map.verticalVideoDistance=loadedData{1}.verticalVideoDistance;

