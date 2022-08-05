function ax = FourChannelInitialization

for i=3:-1:0
ax(4-i).position = [0 i/4 1 1/4];
end

for i = 1:4
       ax(i).label = num2str(i);
end

for i=1:4
    ax(i).status = 2;
end
   

for i=1:4
    ax(i).channel = i - 1;
end


save FourChannelInitialization ax;