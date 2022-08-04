function RegressOutRunning(Fc2,Fs,ybin,forward,rotation,saveDir,MouseID,env,report)
%% Average Over all runs
% dF_ups=Fc2(useP_inds,1);
% for_ups=forward(useP_inds,1);
% rot_ups=rotation(useP_inds,1);
% tot_ups=sqrt(for_ups.^2+rot_ups.^2);
%
% for_reg=[ones(size(for_ups)),for_ups];
% for_out=[ones(size(forward(:,1))),forward(:,1)];
% tot_reg=[ones(size(for_ups)),tot_ups];
% tot_out=[ones(size(forward(:,1))),sqrt(forward(:,1).^2+rotation(:,1).^2)];
%
% slopeFor=for_reg\(dF_ups);
% slopeTot=tot_reg\(dF_ups);
%
% yFor=for_out*slopeFor;
% yTot=tot_out*slopeTot;
%
%
% yForReg=for_reg*slopeFor;
% yTotReg=tot_reg*slopeTot;
%
% figure
% s1=subplot(2,2,1);
% hold on
% plot(Fc2(:,1));
% plot(yTot)
% plot(y/180)
% s2=subplot(2,2,3);
% plot(Fc2(:,1)-yTot);
% s3=subplot(2,2,2);
% plot(dF_ups)
% hold on
% plot(yTotReg)
% plot(y(useP_inds)/180)
% s4=subplot(2,2,4);
% plot(dF_ups-yTotReg)
% linkaxes([s3,s4],'x');
% linkaxes([s1, s2],'x');
% suptitle('Total, Linear Scale')
%
% figure
% s1=subplot(2,2,1);
% hold on
% plot(Fc2(:,1));
% plot(yFor)
% plot(y/180)
% s2=subplot(2,2,3);
% plot(Fc2(:,1)-yFor);
% s3=subplot(2,2,2);
% plot(dF_ups)
% hold on
% plot(yForReg)
% plot(y(useP_inds)/180)
% s4=subplot(2,2,4);
% plot(dF_ups-yForReg)
% linkaxes([s3,s4],'x');
% linkaxes([s1, s2],'x');
% suptitle('Forward, Linear Scale')
% log_for_reg=[ones(size(for_ups)),log(for_ups)];
% log_for_out=[ones(size(forward(:,1))),log(forward(:,1))];
% log_tot_reg=[ones(size(for_ups)),log(tot_ups)];
% log_tot_out=[ones(size(forward(:,1))),log(sqrt(forward(:,1).^2+rotation(:,1).^2))];
%
% log_slopeFor=log_for_reg\(dF_ups);
% log_slopeTot=log_tot_reg\(dF_ups);
%
% log_yFor=log_for_out*log_slopeFor;
% log_yTot=log_tot_out*log_slopeTot;
%
%
% log_yForReg=log_for_reg*log_slopeFor;
% log_yTotReg=log_tot_reg*log_slopeTot;
%
% figure
% s1=subplot(2,2,1);
% hold on
% plot(Fc2(:,1));
% plot(log_yFor)
% plot(y/180)
% s2=subplot(2,2,3);
% plot(Fc2(:,1)-log_yFor);
% s3=subplot(2,2,2);
% plot(dF_ups)
% hold on
% plot(log_yForReg)
% plot(y(useP_inds)/180)
% s4=subplot(2,2,4);
% plot(dF_ups-log_yForReg)
% linkaxes([s3,s4],'x');
% linkaxes([s1, s2],'x');
% suptitle('Forward, Log Scale')
%
% figure
% s1=subplot(2,2,1);
% hold on
% plot(Fc2(:,1));
% plot(log_yTot)
% plot(y/180)
% s2=subplot(2,2,3);
% plot(Fc2(:,1)-log_yTot);
% s3=subplot(2,2,2);
% plot(dF_ups)
% hold on
% plot(log_yTotReg)
% plot(y(useP_inds)/180)
% s4=subplot(2,2,4);
% plot(dF_ups-log_yTotReg)
% linkaxes([s3,s4],'x');
% linkaxes([s1, s2],'x');
% suptitle('Total, Log Scale')
%% For each individual run
% forward=bsxfun(@rdivide, forward,max(forward));
% rotation=bsxfun(@rdivide, rotation,max(rotation));
% for i=1:max(pL)
%     inds=intersect(intersect(find(pL==i), find(y<90)),find(y>10));
%     indout=intersect(intersect(find(pL==i),find(y<170)),find(y>10));
%     dF_ups=Fc2(inds,1);
%     for_ups=forward(inds,1);
%     rot_ups=rotation(inds,1);
%     tot_ups=sqrt(for_ups.^2+rot_ups.^2);
%
%     for_reg=[ones(size(for_ups)),for_ups];
%     for_out=[ones(size(forward(indout,1))),forward(indout,1)];
%     tot_reg=[ones(size(for_ups)),tot_ups];
%     tot_out=[ones(size(forward(indout,1))),sqrt(forward(indout,1).^2+rotation(indout,1).^2)];
%
%     slopeFor=for_reg\(dF_ups);
%     slopeTot=tot_reg\(dF_ups);
%
%     yFor=for_out*slopeFor;
%     yTot=tot_out*slopeTot;
%
%
%     yForReg=for_reg*slopeFor;
%     yTotReg=tot_reg*slopeTot;
%
%     figure
%     s1=subplot(2,2,1);
%     hold on
%     plot(Fc2(indout,1));
%     plot(yTot)
%     plot(y(indout)/180)
%     s2=subplot(2,2,3);
%     plot(Fc2(indout,1)-yTot);
%     s3=subplot(2,2,2);
%     plot(dF_ups)
%     hold on
%     plot(yTotReg)
%     plot(y(inds)/180)
%     s4=subplot(2,2,4);
%     plot(dF_ups-yTotReg)
%     linkaxes([s3,s4],'x');
%     linkaxes([s1, s2],'x');
%     suptitle('Total, Linear Scale')
%
%     figure
%     s1=subplot(2,2,1);
%     hold on
%     plot(Fc2(indout,1));
%     plot(yFor)
%     plot(y(indout)/180)
%     s2=subplot(2,2,3);
%     plot(Fc2(indout,1)-yFor);
%     s3=subplot(2,2,2);
%     plot(dF_ups)
%     hold on
%     plot(yForReg)
%     plot(y(inds)/180)
%     s4=subplot(2,2,4);
%     plot(dF_ups-yForReg)
%     linkaxes([s3,s4],'x');
%     linkaxes([s1, s2],'x');
%     suptitle('Forward, Linear Scale')
%
% end
%
%%
len1=min(Fs*1,100);
len2=min(Fs*4,200);
for cellnum=1:size(ybin,2)
    y=ybin(:,cellnum);
    Fc2=bsxfun(@rdivide, Fc2, max(Fc2));
    [pL,p,nL,n]=split_paths(ybin(:,cellnum),Fs);
    y=y-min(y);
    y=y/max(y)*180;
    runnums=unique(pL);
    runnums=runnums(2:end);
    meanfor=zeros(len1,length(runnums));
    meanrot=zeros(len1,length(runnums));
    meanFc2=zeros(len1,length(runnums));
    meanforout=zeros(len2,length(runnums));
    meanrotout=zeros(len2,length(runnums));
    meanFc2out=zeros(len2,length(runnums));
    for i=1:length(runnums)
        inds=intersect(intersect(find(pL==runnums(i)), find(y<90)),find(y>10));
        tempfor=forward(inds,cellnum);
        temprot=rotation(inds,cellnum);
        tempFc2=Fc2(inds,cellnum);
        
        downratio=length(tempfor)/len1;
        for jj=1:len1
            if (jj*downratio)<len1
                meanfor(jj,i)=mean(tempfor(round(((jj-1)*downratio+1)):round(downratio*jj)));
                meanrot(jj,i)=mean(temprot(round(((jj-1)*downratio+1)):round(downratio*jj)));
                meanFc2(jj,i)=mean(tempFc2(round(((jj-1)*downratio+1)):round(downratio*jj)));
            else
                meanfor(jj,i)=mean(tempfor(round(((jj-1)*downratio+1)):end));
                meanrot(jj,i)=mean(temprot(round(((jj-1)*downratio+1)):end));
                meanFc2(jj,i)=mean(tempFc2(round(((jj-1)*downratio+1)):end));
            end
        end
        
        indout=intersect(intersect(find(pL==runnums(i)),find(y<170)),find(y>10));
        tempforout=forward(indout,cellnum);
        temprotout=rotation(indout,cellnum);
        tempFc2out=Fc2(indout,cellnum);
        downratio=length(tempforout)/len2;
        for jj=1:len2
            if (jj*downratio)<len2
                meanforout(jj,i)=mean(tempforout(round(((jj-1)*downratio+1)):round(downratio*jj)));
                meanrotout(jj,i)=mean(temprotout(round(((jj-1)*downratio+1)):round(downratio*jj)));
                meanFc2out(jj,i)=mean(tempFc2out(round(((jj-1)*downratio+1)):round(downratio*jj)));
            else
                meanforout(jj,i)=mean(tempforout(round(((jj-1)*downratio+1)):end));
                meanrotout(jj,i)=mean(temprotout(round(((jj-1)*downratio+1)):end));
                meanFc2out(jj,i)=mean(tempFc2out(round(((jj-1)*downratio+1)):end));
            end
        end
    end
    
    mfor=nanmean(meanfor,2);
    mrot=nanmean(meanrot,2);
    mtot=sqrt(mfor.^2+mrot.^2);
    mFc2=nanmean(meanFc2,2);
    mforout=nanmean(meanforout,2);
    mrotout=nanmean(meanrotout,2);
    mtotout=sqrt(mforout.^2+mrotout.^2);
    mFc2outup=nanmean(meanFc2out,2);
    
    for_reg=[ones(size(mfor)),mfor];
    for_out=[ones(size(mforout)),mforout];
    tot_reg=[ones(size(mtot)),mtot];
    tot_out=[ones(size(mtotout)),mtotout];
    
    slopeFor=for_reg\(mFc2);
    slopeTot=tot_reg\(mFc2);
    
    yForup=for_out*slopeFor;
    yTotup=tot_out*slopeTot;
    
    
    yForReg=for_reg*slopeFor;
    yTotReg=tot_reg*slopeTot;
    % figure('units','normalized', 'Position', [.01 .05 .49 .87]);
    % subplot(3,1,1)
    % imagesc(meanfor')
    % title('Forward Speed')
    % subplot(3,1,2)
    % imagesc(sqrt(meanfor.^2+meanrot.^2)')
    % title('Total Speed')
    % subplot(3,1,3)
    % imagesc(meanFc2')
    % title('dF/F')
    % suptitle([MouseID, ' Cell ',num2str(cellnum)])
    % saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Up Running Activity', env,'Cell',num2str(cellnum),   '.jpg']);
    % savefig([saveDir,'\',MouseID,'\', MouseID, 'Up Running Activity ',env,'Cell',num2str(cellnum),   '.fig']);
    
    figure('units','normalized', 'Position', [.01 .05 .49 .87]);
    subplot(2,2,1); plot(mFc2outup); hold on; plot(yTotup); plot(abs(mFc2outup-yTotup));
    legend({'Fc2','Estimated Activity','Difference'})
    title('Regression Upwards with Total Speed')
    subplot(2,2,3); plot(mFc2outup); hold on; plot(yForup);plot(abs(mFc2outup-yForup));
    legend({'Fc2','Estimated Activity','Difference'})
    title('Regression Upwards with Forward Speed')
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Regress Up', env, 'Cell',num2str(cellnum),  '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'Regress Up ',env, 'Cell',num2str(cellnum),  '.fig']);
    
    %%
    runnums=unique(nL);
    runnums=runnums(2:end);
    meanfor=zeros(len1,length(runnums));
    meanrot=zeros(len1,length(runnums));
    meanFc2=zeros(len1,length(runnums));
    meanforout=zeros(len2,length(runnums));
    meanrotout=zeros(len2,length(runnums));
    meanFc2out=zeros(len2,length(runnums));
    for i=1:length(runnums)
        inds=intersect(intersect(find(nL==runnums(i)), find(y>90)),find(y<170));
        tempfor=forward(inds,cellnum);
        temprot=rotation(inds,cellnum);
        tempFc2=Fc2(inds,cellnum);
        
        downratio=length(tempfor)/len1;
        for jj=1:len1
            if (jj*downratio)<len1
                meanfor(jj,i)=mean(tempfor(round(((jj-1)*downratio+1)):round(downratio*jj)));
                meanrot(jj,i)=mean(temprot(round(((jj-1)*downratio+1)):round(downratio*jj)));
                meanFc2(jj,i)=mean(tempFc2(round(((jj-1)*downratio+1)):round(downratio*jj)));
            else
                meanfor(jj,i)=mean(tempfor(round(((jj-1)*downratio+1)):end));
                meanrot(jj,i)=mean(temprot(round(((jj-1)*downratio+1)):end));
                meanFc2(jj,i)=mean(tempFc2(round(((jj-1)*downratio+1)):end));
            end
        end
        
        indout=intersect(intersect(find(nL==runnums(i)),find(y<170)),find(y>10));
        tempforout=forward(indout,cellnum);
        temprotout=rotation(indout,cellnum);
        tempFc2out=Fc2(indout,cellnum);
        downratio=length(tempforout)/len2;
        for jj=1:len2
            if (jj*downratio)<len2
                meanforout(jj,i)=mean(tempforout(round(((jj-1)*downratio+1)):round(downratio*jj)));
                meanrotout(jj,i)=mean(temprotout(round(((jj-1)*downratio+1)):round(downratio*jj)));
                meanFc2out(jj,i)=mean(tempFc2out(round(((jj-1)*downratio+1)):round(downratio*jj)));
            else
                meanforout(jj,i)=mean(tempforout(round(((jj-1)*downratio+1)):end));
                meanrotout(jj,i)=mean(temprotout(round(((jj-1)*downratio+1)):end));
                meanFc2out(jj,i)=mean(tempFc2out(round(((jj-1)*downratio+1)):end));
            end
        end
    end
    
    mfor=nanmean(meanfor,2);
    mrot=nanmean(meanrot,2);
    mtot=sqrt(mfor.^2+mrot.^2);
    mFc2=nanmean(meanFc2,2);
    mforout=nanmean(meanforout,2);
    mrotout=nanmean(meanrotout,2);
    mtotout=sqrt(mforout.^2+mrotout.^2);
    mFc2outdown=nanmean(meanFc2out,2);
    
    for_reg=[ones(size(mfor)),mfor];
    for_out=[ones(size(mforout)),mforout];
    tot_reg=[ones(size(mtot)),mtot];
    tot_out=[ones(size(mtotout)),mtotout];
    
    slopeFor=for_reg\(mFc2);
    slopeTot=tot_reg\(mFc2);
    
    yFordown=for_out*slopeFor;
    yTotdown=tot_out*slopeTot;
    
    
    yForReg=for_reg*slopeFor;
    yTotReg=tot_reg*slopeTot;
    % figure('units','normalized', 'Position', [.5 .05 .5 .87]);
    % subplot(3,1,1)
    % imagesc(meanfor')
    % title('Forward Speed')
    % subplot(3,1,2)
    % imagesc(sqrt(meanfor.^2+meanrot.^2)')
    % title('Total Speed')
    % subplot(3,1,3)
    % imagesc(meanFc2')
    % title('dF/F')
    % suptitle([MouseID, ' Cell ',num2str(cellnum)])
    % saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Down Running Activity', env,  'Cell',num2str(cellnum), '.jpg']);
    % savefig([saveDir,'\',MouseID,'\', MouseID, 'Down Running Activity ',env,  'Cell',num2str(cellnum), '.fig']);
    
    
    % figure('units','normalized', 'Position', [.5 .05 .5 .87]);
    
    subplot(2,2,2); plot(mFc2outdown); hold on; plot(yTotdown); plot(abs(mFc2outdown-yTotdown));
    legend({'Fc2','Estimated Activity','Difference'})
    title('Regression Downwards with Total Speed')
    subplot(2,2,4); plot(mFc2outdown); hold on; plot(yFordown);plot(abs(mFc2outdown-yFordown));
    legend({'Fc2','Estimated Activity','Difference'})
    title('Regression Downwards with Forward Speed')
    import mlreportgen.dom.*;
    suptitle([MouseID, ' Cell ',num2str(cellnum)])
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Regress ', env,   'Cell',num2str(cellnum),'.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'Regress ',env,   'Cell',num2str(cellnum),'.fig']);
    
    
end

if exist('report','var')
    import mlreportgen.dom.*;
    regTable=Table;
    regRow=TableRow;
    for cellnum=1:size(ybin,2)
        temp=Image([saveDir,'\',MouseID,'\', MouseID, 'Regress ', env,   'Cell',num2str(cellnum),'.jpg']);
        append(regRow,TableEntry(temp));
    end
    append(regTable,regRow);
    append(report,regTable);
end