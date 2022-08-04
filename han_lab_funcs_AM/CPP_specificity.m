% uregionPref{1}=zeros(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],2),3);
% uregionPref{1}(:,1)=mean([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}])';
% uregionPref{1}(:,2)=mean([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])';
% uregionPref{1}(:,3)=mean([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}])';
% semregionPref{1}=zeros(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],2),3);
% semregionPref{1}(:,1)=std([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],1))';
% semregionPref{1}(:,2)=std([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}],1))';
% semregionPref{1}(:,3)=std([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}],1))';
% 
% uregionPref{2}=zeros(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],2),3);
% uregionPref{2}(:,1)=mean([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}])';
% uregionPref{2}(:,2)=mean([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])';
% uregionPref{2}(:,3)=mean([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}])';
% semregionPref{2}=zeros(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],2),3);
% semregionPref{2}(:,1)=std([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],1))';
% semregionPref{2}(:,2)=std([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}],1))';
% semregionPref{2}(:,3)=std([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}],1))';


%Positive y is b, negative y is A
maxh=876.2818;
m=2;
curr=1;
useinds=1:round(5.2*60*15);
yb=(((mouseCPP(m).ybinned{curr}{1}(useinds,1))/maxh)+1)*1.5;
down=yb<1;
up=yb>2;
regionTime{3}{1}=sum(up)/(sum(up)+sum(down));
regionTime{3}{2}=sum(down)/(sum(up)+sum(down));
uregionPref{3}=zeros(size([mouseCPP(m).FPyrs{curr}{1}],2),3);
uregionPref{3}(:,1)=mean(mouseCPP(m).FPyrs{curr}{1}(up,:))';
uregionPref{3}(:,2)=mean([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])';
uregionPref{3}(:,3)=mean(mouseCPP(m).FPyrs{curr}{1}(down,:))';
semregionPref{3}=zeros(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],2),3);
% semregionPref{3}(:,1)=std([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],1))';
% semregionPref{3}(:,2)=std([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}],1))';
% semregionPref{3}(:,3)=std([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}],1))';

m=2;
curr=2;
yb=(((mouseCPP(m).ybinned{curr}{1}(useinds,1))/maxh)+1)*1.5;
down=yb<1;
up=yb>2;
regionTime{4}{1}=sum(up)/(sum(up)+sum(down));
regionTime{4}{2}=sum(down)/(sum(up)+sum(down));
uregionPref{4}=zeros(size([mouseCPP(m).FPyrs{curr}{1}],2),3);
uregionPref{4}(:,1)=mean(mouseCPP(m).FPyrs{curr}{1}(up,:))';
uregionPref{4}(:,2)=mean([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])';
uregionPref{4}(:,3)=mean(mouseCPP(m).FPyrs{curr}{1}(down,:))';
semregionPref{4}=zeros(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],2),3);
% semregionPref{4}(:,1)=std([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],1))';
% semregionPref{4}(:,2)=std([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}],1))';
% semregionPref{4}(:,3)=std([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}],1))';


% uregionPref{5}=zeros(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],2),3);
% uregionPref{5}(:,1)=mean([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}])';
% uregionPref{5}(:,2)=mean([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])';
% uregionPref{5}(:,3)=mean([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}])';
% semregionPref{5}=zeros(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],2),3);
% semregionPref{5}(:,1)=std([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],1))';
% semregionPref{5}(:,2)=std([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}],1))';
% semregionPref{5}(:,3)=std([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}],1))';
% 
% uregionPref{6}=zeros(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],2),3);
% uregionPref{6}(:,1)=mean([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}])';
% uregionPref{6}(:,2)=mean([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])';
% uregionPref{6}(:,3)=mean([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}])';
% semregionPref{6}=zeros(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],2),3);
% semregionPref{6}(:,1)=std([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{1} mouseCPP(m).IntActivity{curr}{1}],1))';
% semregionPref{6}(:,2)=std([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{2} mouseCPP(m).IntActivity{curr}{2}],1))';
% semregionPref{6}(:,3)=std([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}])/sqrt(size([mouseCPP(m).PyrActivity{curr}{3} mouseCPP(m).IntActivity{curr}{3}],1))';

uregionPref{3}(uregionPref{3}<0)=0;
uregionPref{4}(uregionPref{4}<0)=0;

regionSpec{6}=0;
% regionSpec{1}=(uregionPref{1}(:,1)-uregionPref{1}(:,3))./(uregionPref{1}(:,1)+uregionPref{1}(:,3));
% regionSpec{2}=(uregionPref{2}(:,1)-uregionPref{2}(:,3))./(uregionPref{2}(:,1)+uregionPref{2}(:,3));
regionSpec{3}=(uregionPref{3}(:,1)-uregionPref{3}(:,3))./(uregionPref{3}(:,1)+uregionPref{3}(:,3));
regionSpec{4}=(uregionPref{4}(:,1)-uregionPref{4}(:,3))./(uregionPref{4}(:,1)+uregionPref{4}(:,3));
% regionSpec{5}=(uregionPref{5}(:,1)-uregionPref{5}(:,3))./(uregionPref{5}(:,1)+uregionPref{5}(:,3));
% regionSpec{6}=(uregionPref{6}(:,1)-uregionPref{6}(:,3))./(uregionPref{6}(:,1)+uregionPref{6}(:,3));

% figure; imagesc(uregionPref{1}')
% figure; imagesc(bsxfun(@rdivide,uregionPref{1},max(uregionPref{1},[],2))')
% figure; histogram(regionSpec{1})
% figure; ksdensity(regionSpec{1})
% figure; stem(sort(regionSpec{1}))
% 
% 
% figure; imagesc(uregionPref{2}')
% figure; imagesc(bsxfun(@rdivide,uregionPref{2},max(uregionPref{2},[],2))')
% figure; histogram(regionSpec{2})
% figure; ksdensity(regionSpec{2})
% figure; stem(sort(regionSpec{2}))



figure;bar([regionTime{3}{1} regionTime{3}{2}])
figure; imagesc(uregionPref{3}')
figure; imagesc(bsxfun(@rdivide,uregionPref{3},max(uregionPref{3},[],2))')
figure; histogram(regionSpec{3})
figure; ksdensity(regionSpec{3})
stPos=sort(regionSpec{3});
stNeg=sort(regionSpec{3});
stPos(stPos<0)=nan;
stNeg(stNeg>0)=nan;
figure; stem(stNeg,'color','red');
hold on; stem(stPos,'color','b');
figure; stem(sort(regionSpec{3}))
bar([(sum(regionSpec{3}>0)/length(regionSpec{3})),(sum(regionSpec{3}<0)/length(regionSpec{3}))])


figure;bar([regionTime{4}{1} regionTime{4}{2}])
figure; imagesc(uregionPref{4}')
figure; imagesc(bsxfun(@rdivide,uregionPref{4},max(uregionPref{4},[],2))')
figure; histogram(regionSpec{4})
figure; ksdensity(regionSpec{4})
stPos=sort(regionSpec{4});
stNeg=sort(regionSpec{4});
stPos(stPos<0)=nan;
stNeg(stNeg>0)=nan;
figure; stem(stNeg,'color','red');
hold on; stem(stPos,'color','b');
figure; stem(sort(regionSpec{4}))
figure; bar([(sum(regionSpec{4}>0)/length(regionSpec{4})),(sum(regionSpec{4}<0)/length(regionSpec{4}))])

% figure; imagesc(uregionPref{5}')
% figure; imagesc(bsxfun(@rdivide,uregionPref{5},max(uregionPref{5},[],2))')
% figure; histogram(regionSpec{5})
% figure; ksdensity(regionSpec{5})
% figure; stem(sort(regionSpec{5}))
% 
% figure; imagesc(uregionPref{6}')
% figure; imagesc(bsxfun(@rdivide,uregionPref{6},max(uregionPref{6},[],2))')
% figure; histogram(regionSpec{6})
% figure; ksdensity(regionSpec{6})
% figure; stem(sort(regionSpec{6}))