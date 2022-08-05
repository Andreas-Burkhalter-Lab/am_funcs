function ax = BundleInitialization

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

for i = 0:63
       ax(i+1).label = num2str(i);
end
%ax(1).label = 'Stim';



for i=1:64
    ax(i).status = 2;
end

for i=1:64
    ax(i).channel = i - 1;
end

save BundleInitialization ax   