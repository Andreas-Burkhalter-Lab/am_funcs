function plot_xcorr_lags(mouse,varargin)
if nargin>1
    mouseNums=varargin{1};
    expnum=varargin{2};
    Fs=varargin{3};
    if nargin>4
        saveDir=[varargin{4},'\Lags'];
    else     saveDir='F:\MA Data\Interneurons\PubFigs\Lags\';
    end
    
else
    mouseNums=1:length(mouse);
    for m=mouseNums
        expnum{m}=1:length(mouse(m).Falls);
    end
    Fs=4*ones(size(mouseNums));
    saveDir='F:\MA Data\Interneurons\PubFigs\Lags\';
    
end
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
forcorrs=[];
totcorrs=[];
vrcorrs=[];
prta=[];
%prtas shape: cells 7*3
%prta{1}=6x8
for m=mouseNums
    for e=expnum{m}
        m
        e
        f=(mouse(m).mfballout(e));
        t=mouse(m).mtballout(e);
        v=mouse(m).mvrout(e);
        p=mouse(m).mPRTAout(e);
        forcorrs=[forcorrs,f{1}];
        totcorrs=[totcorrs,t{1}];
        vrcorrs=[vrcorrs,v{1}];
        prta=[prta,p{1}];
    end
end

plot_instance(forcorrs,'Forward');
plot_instance(totcorrs,'Total Ball');
plot_instance(vrcorrs,'VR');
plot_prta_instance(prta);

    function plot_instance(info,name)
        figure
        info(2,isnan(info(1,:)))=nan;
        histogram(info(2,:),50)
        title([name, ' lag of max xcorr'])
        saveas(gca,[saveDir,name,'.jpg']);
        savefig([saveDir,name,'.fig']);
    end
    function plot_prta_instance(info)
        %Max and mins of PRTAs
        %To use the 8s PRTA information its 8+(n-1)*21
        temp=[];
        for ii=1:(size(info,2)/3)
            temp=[temp,info{8+(ii-1)*21}];
            numcells(ii)=size(info{8+(ii-1)*21},2);
        end;
        figure
        prtapk1=min(temp(2,:),temp(4,:));
        prtapk2=max(temp(2,:),temp(4,:));
        prtapk3=temp(6,:);
        subplot(3,1,1)
        histogram(prtapk1,50)
        title('First Peak Time in PRTA traces')
        subplot(3,1,2)
        histogram(prtapk2,50)
        title('Seond Peak Time in PRTA traces')
        subplot(3,1,3)
        histogram(prtapk3,50)
        title('Minimum Time in PRTA traces')
        saveas(gca,[saveDir,'PRTA max and min times','.jpg']);
        savefig([saveDir,'PRTA max and min times','.fig']);
        %Max and mins of forward running PRTAs
        temp=[];
        for ii=1:(size(info,2)/3)
            temp=[temp,(info{9+(ii-1)*21})'];
        end;
        figure
        forpk1=min(temp(2,:),temp(4,:));
        forpk2=max(temp(2,:),temp(4,:));
        forpk3=temp(6,:);
        subplot(3,1,1)
        histogram(forpk1,50)
        title('First Peak Time in forward running PRTA traces')
        subplot(3,1,2)
        histogram(forpk2,50)
        title('Seond Peak Time in forward running PRTA traces')
        subplot(3,1,3)
        histogram(forpk3,50)
        title('Minimum Time in forward running PRTA traces')
        saveas(gca,[saveDir,'forward PRTA max and min times','.jpg']);
        savefig([saveDir,'forward PRTA max and min times','.fig']);
        %Max and mins of total running PRTAs
        temp=[];
        for ii=1:(size(info,2)/3)
            temp=[temp,(info{10+(ii-1)*21})'];
        end;
        figure
        totpk1=min(temp(2,:),temp(4,:));
        totpk2=max(temp(2,:),temp(4,:));
        totpk3=temp(6,:);
        subplot(3,1,1)
        histogram(totpk1,50)
        title('First Peak Time in total running PRTA traces')
        subplot(3,1,2)
        histogram(totpk2,50)
        title('Seond Peak Time in running PRTA traces')
        subplot(3,1,3)
        histogram(totpk3,50)
        title('Minimum Time in running PRTA traces')
        saveas(gca,[saveDir,'total PRTA max and min times','.jpg']);
        savefig([saveDir,'total PRTA max and min times','.fig']);
        
        %Max lag times of forward running dF PRTA xcorr
        tempprta=[];
        for ii=1:(size(info,2)/3)
            tempprta=[tempprta,(info{12+(ii-1)*21})];
        end
        figure
        pk1=tempprta(2,:);
        subplot(2,1,1)
        histogram(pk1,50)
        title('Forward running dF PRTA xcorr lag time')
        temp=[];
        for ii=1:(size(info,2)/3)
            temp=[temp,(info{13+(ii-1)*21})];
        end
        pk1=temp(2,:);
        subplot(2,1,2)
        histogram(pk1,50)
        title('Total running dF PRTA xcorr lag time')
        
        saveas(gca,[saveDir,'running PRTA xcorr lag times','.jpg']);
        savefig([saveDir,'running PRTA xcorr lag times','.fig']);
        
%         %Difference between forward running peaks and dF peaks PRTA
%         %To use the 8s PRTA information its 8+(n-1)*21
%         tempprta=[];
%         for ii=1:(size(info,2)/3)
%             temp=info{8+(ii-1)*21};
%             tempfor=(info{9+(ii-1)*21});
%             temp=bsxfun(@minus,temp,tempfor);
%             tempprta=[tempprta,temp];
%         end;
%         figure
%         prtapk1=min(tempprta(2,:),tempprta(4,:));
%         prtapk2=max(tempprta(2,:),tempprta(4,:));
%         prtapk3=tempprta(6,:);
%         subplot(3,1,1)
%         histogram(prtapk1,50)
%         title('Difference in first peak PRTA, forward')
%         subplot(3,1,2)
%         histogram(prtapk2,50)
%         title('Difference in Second peak PRTA, forward')
%         subplot(3,1,3)
%         histogram(prtapk3,50)
%         title('Difference in Third peak PRTA, forward')
%         saveas(gca,[saveDir,'difference between PRTA and forward peaks','.jpg']);
%         savefig([saveDir,'difference between PRTA and forward peaks','.fig']);
%         
%                 %Difference between total running peaks and dF peaks PRTA
%         %To use the 8s PRTA information its 8+(n-1)*21
%         tempprta=[];
%         for ii=1:(size(info,2)/3)
%             temp=info{8+(ii-1)*21};
%             tempfor=(info{10+(ii-1)*21})';
%             temp=bsxfun(@minus,temp,tempfor);
%             tempprta=[tempprta,temp];
%         end;
%         figure
%         prtapk1=min(tempprta(2,:),tempprta(4,:));
%         prtapk2=max(tempprta(2,:),tempprta(4,:));
%         prtapk3=tempprta(6,:);
%         subplot(3,1,1)
%         histogram(prtapk1,50)
%         title('Difference in first peak PRTA, total')
%         subplot(3,1,2)
%         histogram(prtapk2,50)
%         title('Difference in Second peak PRTA, total')
%         subplot(3,1,3)
%         histogram(prtapk3,50)
%         title('Difference in Third peak PRTA, total')
%         saveas(gca,[saveDir,'difference between PRTA and total peaks','.jpg']);
%         savefig([saveDir,'difference between PRTA and total peaks','.fig']);

        

    end
end