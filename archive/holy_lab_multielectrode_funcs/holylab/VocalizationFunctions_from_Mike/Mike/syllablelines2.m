function [linesOut,numLines,pfOut] = syllablelines2(snip,thresh,scale,f)
%% Make an image

x = full(abs(snip));
x(x==0) = thresh;
x = imfilter(x,fspecial('gaussian',[10,10],2),'replicate');
x = 10.*log10(x.^2);
x = roundn(x,-1);
x = x./min(x(:));

%% Label components

z = bwlabeln(x>min(x(:)));

%% Look at each component

components = unique(z(z>0));

for i = 1:numel(components)
    z2 = x;
    k = find(z~=components(i));
    z2(k) = 1;
    [maxVals,linesOut.line{i}] = max(z2);
    linesOut.line{i} = f(linesOut.line{i});
    linesOut.line{i}(maxVals <= scale) = 0;
    linesOut.mf(i) = mean(linesOut.line{i}(linesOut.line{i}>0));
end

k = find(isnan(linesOut.mf));
linesOut.line(k) = [];
linesOut.mf(k) = [];

for i = 1:numel(linesOut.mf)
    for j = 1:numel(linesOut.mf)
        linesOut.overlap(i,j) = numel(intersect(find(linesOut.line{i}>0),find(linesOut.line{j}>0)));
    end
end


numLines = numel(linesOut.line);

%% Calculate pf

[maxVals,pfOut] = max(x);
pfOut = f(pfOut);
pfOut(maxVals <= scale) = 0;
end