function RFfitting2case
%% case 1
data=xlsread('BookJI.xlsx','Sheet1');
data(any(isnan(data),2),:)=[];
[~, index2] = sort(data(:,1));
newdata = data(index2,:);

scatter(newdata(:,1),newdata(:,2),'r.')
xlabel('Cell to Cell Distance (\mum)')
ylabel('Y ARF Overlap (degree)')
title('Receptive Field Sigmoid Fitting')

% Sigmoidal Fitting
x=newdata(:,1);
y=newdata(:,2);
f = @(p,x) p(1) + p(2) ./ (1+exp(p(3)*(x-p(4))));
p = nlinfit(x,y,f,[-40 40 .001 200])
fprintf('Sigmoidal Fitting Equation=\n %2.2f + |       %2.2f       | \n          |-------------------|\n          |1+e^%.4f(x-%3.1f)|\n',p(1),p(2),p(3),p(4));
hold on
t=-200:1000;
z=p(1)+p(2)./(1+exp(p(3)*(t'-p(4))));
plot(t,z,'r')
baseline=zeros(1,length(t));
plot(t,baseline,'k')
xpt=find(z<=.01,1,'first');
fprintf('Zero_Crossing(y=0) at x=%3.0f um \n',t(xpt))
y0=0;
plot(t(xpt),y0,'r^')
% std
leng=length(newdata)-1;
brk=zeros(1,leng);
for n=1:leng
    brk(n)=newdata(n+1,1)-newdata(n,1);
end

pt = find(brk~=0);
lmatrix=length(pt)+1;

stdmatrix=zeros(1,lmatrix);
meanmatrix=zeros(1,lmatrix);
sematrix=zeros(1,lmatrix);
for n=1:lmatrix
    if n==1
        stdmatrix(n)=std(newdata(1:pt(n),2));
        sematrix(n)=std(newdata(1:pt(n),2))/sqrt(length(newdata(1:pt(n),2)));
        meanmatrix(n)=mean(newdata(1:pt(n),2));
    elseif n==lmatrix
        stdmatrix(n)=std(newdata(pt(n-1)+1:end,2));
        meanmatrix(n)=mean(newdata(pt(n-1)+1:end,2));
        sematrix(n)=std(newdata(pt(n-1)+1:end,2))/sqrt(length(newdata(pt(n-1)+1:end,2)));
    else
        stdmatrix(n)=std(newdata(pt(n-1)+1:pt(n),2));
        meanmatrix(n)=mean(newdata(pt(n-1)+1:pt(n),2));        
        sematrix(n)=std(newdata(pt(n-1):pt(n),2))/sqrt(length(newdata(pt(n-1):pt(n),2)));
    end
end
Group1std=stdmatrix
Group1se=sematrix
for n=1:lmatrix
    if n==lmatrix
        X(n)=newdata(end,1);
    else
        X(n)=newdata(pt(n),1);
    end
end
Y=p(1)+p(2)./(1+exp(p(3)*(X-p(4))));

bot=zeros(lmatrix,5);
top=zeros(lmatrix,5);
line=zeros(lmatrix,5);
dat=zeros(lmatrix,5);
mid=zeros(lmatrix,5);
for n=1:lmatrix
    for m=1:5
        mid(n,m)=Y(n)-sematrix(n)+sematrix(n)*(m-1)/2;
    end
    dat(n,:)=X(n);
    plot(dat(n,:),mid(n,:),'r','linewidth',1.5)
end
for n=1:lmatrix
    for m=1:5
        line(n,m)=X(n)-10+20*(m-1)/4;
    end
    bot(n,:)=Y(n)-sematrix(n);
    top(n,:)=Y(n)+sematrix(n);
    plot(line(n,:),bot(n,:),'r','linewidth',1.5)
    plot(line(n,:),top(n,:),'r','linewidth',1.5)    
end

%% case 2
data2=xlsread('BookJI2.xlsx','Sheet1');
data2(any(isnan(data2),2),:)=[];
[~, index3] = sort(data2(:,1));
newdata2 = data2(index3,:);

scatter(newdata2(:,1),newdata2(:,2),'g+')
xlabel('Cell to Cell Distance (\mum)')
ylabel('Y ARF Overlap (degree)')
title('Receptive Field Sigmoid Fitting')

% Sigmoidal Fitting
x2=newdata2(:,1);
y3=newdata2(:,2);
f2 = @(p,x2) p(1) + p(2) ./ (1+exp(p(3)*(x2-p(4))));
p2 = nlinfit(x2,y3,f2,[-40 40 .001 200])
fprintf('Sigmoidal Fitting Equation=\n %2.2f + |       %2.2f       | \n          |-------------------|\n          |1+e^%.4f(x-%3.1f)|\n',p(1),p(2),p(3),p(4));
hold on
t=-200:1000;
z2=p2(1)+p2(2)./(1+exp(p2(3)*(t'-p2(4))));
plot(t,z2,'g')
xpt2=find(z2<=.01,1,'first');
fprintf('Zero_Crossing(y=0) at x=%3.0f um \n',t(xpt2))
plot(t(xpt2),y0,'g^')
% std
leng2=length(newdata2)-1;
brk2=zeros(1,leng2);
for n=1:leng2
    brk2(n)=newdata2(n+1,1)-newdata2(n,1);
end

pt2 = find(brk2~=0);
lmatrix2=length(pt2)+1;

stdmatrix2=zeros(1,lmatrix2);
meanmatrix2=zeros(1,lmatrix2);
sematrix2=zeros(1,lmatrix2);
for n=1:lmatrix2
    if n==1
        stdmatrix2(n)=std(newdata2(1:pt2(n),2));
        sematrix2(n)=std(newdata2(1:pt2(n),2))/sqrt(length(newdata2(1:pt2(n),2)));
        meanmatrix2(n)=mean(newdata2(1:pt2(n),2));
    elseif n==lmatrix2
        stdmatrix2(n)=std(newdata2(pt2(n-1)+1:end,2));
        meanmatrix2(n)=mean(newdata2(pt2(n-1)+1:end,2));
        sematrix2(n)=std2(newdata2(pt2(n-1)+1:end,2))/sqrt(length(newdata2(pt2(n-1)+1:end,2)));
    else
        stdmatrix2(n)=std(newdata2(pt2(n-1)+1:pt2(n),2));
        meanmatrix2(n)=mean(newdata2(pt2(n-1)+1:pt2(n),2));
        sematrix2(n)=std2(newdata2(pt2(n-1):pt2(n),2))/sqrt(length(newdata2(pt2(n-1):pt2(n),2)));
    end
end
Group2std=stdmatrix2
Group2se=sematrix2
X2=zeros(1,lmatrix2);
for n=1:lmatrix2
    if n==lmatrix2
        X2(n)=newdata2(end,1);
    else
        X2(n)=newdata2(pt2(n),1);
    end
end
Y2=p2(1)+p2(2)./(1+exp(p2(3)*(X2-p2(4))));

bot2=zeros(lmatrix2,5);
top2=zeros(lmatrix2,5);
line2=zeros(lmatrix2,5);
dat2=zeros(lmatrix2,5);
mid2=zeros(lmatrix2,5);
for n=1:lmatrix2
    for m=1:5
        mid2(n,m)=Y2(n)-sematrix2(n)+sematrix2(n)*(m-1)/2;
    end
    dat2(n,:)=X2(n);
    plot(dat2(n,:),mid2(n,:),'g','linewidth',1.5)
end
for n=1:lmatrix2
    for m=1:5
        line2(n,m)=X2(n)-10+20*(m-1)/4;
    end
    bot2(n,:)=Y2(n)-sematrix2(n);
    top2(n,:)=Y2(n)+sematrix2(n);
    plot(line2(n,:),bot2(n,:),'g','linewidth',1.5)
    plot(line2(n,:),top2(n,:),'g','linewidth',1.5)    
end

%% +/- 1 STD
[indexx]=find(abs(newdata(:,2))<=1);
xmatrix=zeros(length(indexx),1);
for o=1:length(indexx)
   xmatrix(o,1)=newdata(indexx(o),1);
end
xstd=std(xmatrix);
fprintf('+/- 1 Value Range STD = %3.3f\n', xstd)
[indexx]=find(abs(newdata2(:,2))<=1);
xmatrix=zeros(length(indexx),1);
for o=1:length(indexx)
   xmatrix(o,1)=newdata2(indexx(o),1);
end
xstd=std(xmatrix);
fprintf('+/- 1 Value Range STD = %3.3f\n', xstd)

% Calculating R^2 value and AnovaP
y2=p(1)+p(2)./(1+exp(p(3)*(x-p(4))));
R2=rsquare(y,y2);
fprintf('R-Sqaure Value = %1.4f\n', R2)
p=anova1(y,y2);
y4=p2(1)+p2(2)./(1+exp(p2(3)*(x2-p2(4))));
R2=rsquare(y3,y4);
fprintf('R-Sqaure Value = %1.4f\n', R2)
p=anova1(y3,y4);


%nonparametric Mann-Whitney U test / parametric
[p,h,stats] = ranksum(z,z2);
STATS=mwwtest(z,z2);

fprintf('Nonparametric Mann-Whitney/Wilcoxon Test for 2 curve fittings \n Mann-Whitney U = %8.3f \n WilcoxonW = %10.3f \n Z = %1.3f \n Asymp.Sig.(2-tailed) = %1.3f \n',STATS.U,STATS.T,stats.zval,p)

% [p,h,stats] = signrank(z,z2)
[p,h,stats] = ranksum(y,y3);
STATS=mwwtest(y,y3);

fprintf('Nonparametric Mann-Whitney/Wilcoxon Test for 2 exact group data\n Mann-Whitney U = %8.3f \n WilcoxonW = %10.3f \n Z = %1.3f \n Asymp.Sig.(2-tailed) = %1.3f \n',STATS.U,STATS.T,stats.zval,p)
end
