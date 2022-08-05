function [syllObj,objMf,numObj,objOverlap,syllPf,syllPfadj,x] = syllableObjects(snip,f,minHeight)
%% Make an image

x = full(abs(snip));
x2 = imfilter(x,fspecial('gaussian',[25 25],5),'replicate');


%% Label components
if max(max(x2)) < minHeight
    z = bwlabeln(x2>0);
else
    z = bwlabeln(x2>minHeight);
end

%% Look at each component
components = unique(z(z>0));
syllObj = cell(1,numel(components));
objMf = zeros(1,numel(components));
for i = 1:numel(components)
    z2 = double(z==components(i));
    [m,n] = find(z2);
    z2(m,n) = x2(m,n);
    [~,syllObj{i}] = max(z2);
    syllObj{i} = f(syllObj{i});
    objMf(i) = mean(syllObj{i}(syllObj{i}>0));
end
k = find(isnan(objMf));
objMf(k) = [];
syllObj(k) = [];
objOverlap = zeros(numel(syllObj),numel(syllObj));
for i = 1:numel(syllObj)
    for j = 1:numel(syllObj)
        objOverlap(i,j) = numel(intersect(find(syllObj{i}>0),find(syllObj{j}>0)));
    end
end
numObj = numel(syllObj);

%% Calculate pf

[~,syllPf] = max(x);
syllPf = f(syllPf);
syllPf(max(x) <= minHeight) = 0; 

%% Use neighbors to fill in zeros for syllPfadj

syllPfadj = syllPf;
k = find(syllPfadj == 0);  
m = find(syllPfadj);
if ~isempty(k) && ~isempty(m)
 for t = 1:numel(k)
    dist = abs(m - k(t));                                    % Subtract index of value-containing entries from the zero entry
    dist = m(find(dist == min(dist(dist>0)),1));             % Index of value containing entries closest. This is usual 1 away but could be more in the case of multiple zeros.
            syllPfadj(k(t)) = syllPfadj(dist);   
 end
elseif isempty(m);
 syllPfadj = NaN(length(syllPf));
end


end