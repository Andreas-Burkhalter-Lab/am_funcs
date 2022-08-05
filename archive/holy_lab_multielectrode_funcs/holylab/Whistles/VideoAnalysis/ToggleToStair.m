function st = ToggleToStair(toggle,ctime,nbehavs)
st = cell(1,nbehavs);
for i = 1:nbehavs
        st{i} = zeros(2,1);
end
ntog = size(toggle,2);
for i = 1:ntog
        indx = toggle(1,i);
        st{indx}(:,end+1) = [~st{indx}(1,end); toggle(2,i)];
end
for i = 1:nbehavs
        st{i}(:,end+1) = [st{i}(1,end); ctime];
end
