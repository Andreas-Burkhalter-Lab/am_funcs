function ax = MEAInitialization


for i=7:-1:0
        ax(8-i).position = [0 i/8 1/8 1/8];
        ax(16-i).position = [1/8 i/8 1/8 1/8];
        ax(24-i).position = [2/8 i/8 1/8 1/8];
        ax(32-i).position = [3/8 i/8 1/8 1/8];
        ax(40-i).position = [4/8 i/8 1/8 1/8];
        ax(48-i).position = [5/8 i/8 1/8 1/8];
        ax(56-i).position = [6/8 i/8 1/8 1/8];
        ax(64-i).position = [7/8 i/8 1/8 1/8];
end


ax(1).label = 'Stim';
for i=12:17
    ax(i-10).label = num2str(i);
end
ax(8).label = 'label2'
for i=21:28
    ax(i-12).label = num2str(i);
end
for i=31:38
    ax(i-14).label = num2str(i);
end
for i=41:48
    ax(i-16).label = num2str(i);
end
for i=51:58
    ax(i-18).label = num2str(i);
end
for i=61:68
    ax(i-20).label = num2str(i);
end
for i=71:78
    ax(i-22).label = num2str(i);
end
ax(57).label = 'label3';
for i=82:87
    ax(i-24).label = num2str(i);
end
ax(64).label = 'label4';

for i=1:64
    ax(i).status = 2;
end
   

for i=1:64
    ax(i).channel = i - 1;
end

save MEAInitialization ax;