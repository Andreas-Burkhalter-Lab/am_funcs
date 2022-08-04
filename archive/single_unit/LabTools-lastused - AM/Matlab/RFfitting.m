function RFfitting
% DoHyun Kim's Summer Project #3, Fitting Script
% Logic development
% 1. Input 3 dimentional data//2 dimentional data
%     x-axis(only one data): cell to cell distance
%     y-axis(more than one data): Y ARF Overlap
%     ,, Need to make a input script
%     ,, NaN Problems
%     
         % M(any(isnan(M),2),:)=[]
         % M(isfinite(M(:, 1)), :)
% 2. Make a Scatter Plot
% 
% 3. Make a fitting (try Sigmoidal // Exponential)
% 
% 4. Statistical Analysis
%     ex) Critical Points, AnovaP, R^2
data=xlsread('BookJI.xlsx','Sheet1');
data(any(isnan(data),2),:)=[];
[~, index2] = sort(data(:,1));
newdata = data(index2,:);

scatter(newdata(:,1),newdata(:,2),'.')
xlabel('Cell to Cell Distance (\mum)')
ylabel('Y ARF Overlap (degree)')
title('Receptive Field Sigmoid Fitting')

% Sigmoidal Fitting
x=newdata(:,1);
y=newdata(:,2);
% P1= min(newdata(:,2)); % lower plateau
% P2= max(newdata(:,2))-min(newdata(:,2)); % Range
% P3= 1; %curvature coefficient
% P4= 0; % BP50.

% y=P1+P2/(1+exp(P3*(x.*1-P4)));

f = @(p,x) p(1) + p(2) ./ (1+exp(p(3)*(x-p(4))));
p = nlinfit(x,y,f,[-40 40 .001 200])
% line(x,f(p,x),'color','r')
fprintf('Sigmoidal Fitting Equation=\n %2.2f + |       %2.2f       | \n          |-------------------|\n          |1+e^%.4f(x-%3.1f)|\n',p(1),p(2),p(3),p(4));

hold on
t=-200:1000;
z=p(1)+p(2)./(1+exp(p(3)*(t'-p(4))));
plot(t,z,'g')
baseline=zeros(1,length(t));
plot(t,baseline,'k')
xpt= find(z<=.01,1,'first');
fprintf('Zero_Crossing(y=0) at x=%3.0f um \n',t(xpt))
y0=0;
plot(t(xpt),y0,'m^')
%% std
leng=length(newdata)-1;
brk=zeros(1,leng);
for n=1:leng
    brk(n)=newdata(n+1,1)-newdata(n,1);
end

pt = find(brk~=0);
lmatrix=length(pt)+1;

stdmatrix=zeros(1,lmatrix);
sematrix=zeros(1,lmatrix);
meanmatrix=zeros(1,lmatrix);
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
        sematrix(n)=std(newdata(pt(n-1)+1:pt(n),2))/sqrt(length(newdata(pt(n-1)+1:pt(n),2)));
    end
end
StandardDeviation=stdmatrix
StandardError=sematrix
X=zeros(1,lmatrix);
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

%% +/- 1 STD
[indexx]=find(abs(newdata(:,2))<=1);
xmatrix=zeros(length(indexx),1);
for o=1:length(indexx)
   xmatrix(o,1)=newdata(indexx(o),1);
end
xstd=std(xmatrix);
fprintf('+/- 1 Value Range STD = %3.3f\n', xstd)

%% Calculating R^2 value and AnovaP
y2=p(1)+p(2)./(1+exp(p(3)*(x-p(4))));
R2=rsquare(y,y2);
fprintf('R-Sqaure Value = %1.4f\n', R2)
p=anova1(y,y2);
end

