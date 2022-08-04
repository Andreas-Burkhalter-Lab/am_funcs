function interneuron_frequencies(Fall,Fs, rewards, pos, speedv,speedy,varargin)
if nargin>6
    saveDir=varargin{1};
    MouseID=varargin{2};
    if nargin>8
        env=varargin{3};
    else env='';
    end
else
saveDir='F:\MA Data\Interneurons\';
MouseID='This Mouse';
end
pos=bsxfun(@minus,pos,min(pos,[],1));
lsig=zeros(size(Fall));
hsig=lsig;
gsig=hsig;
asig=gsig;
betsig=asig;
thesig=betsig;
delsig=thesig;
num_cells=size(Fall,2);
times=linspace(1,length(Fall)/Fs,length(Fall));

for cell=1:num_cells
    lp=designfilt('lowpassfir','FilterOrder',128, ...
        'PassbandFrequency',.2,'StopbandFrequency',.4,...
        'SampleRate',Fs);
    lsig(:,cell)=filtfilt(lp,Fall(:,cell));
    del=mean(grpdelay(lp));
    hp=designfilt('highpassfir','FilterOrder',128, ...
        'PassbandFrequency',.2,'StopbandFrequency',.1,...
        'SampleRate',Fs);
    if Fs>=10
        delp=designfilt('bandpassfir','FilterOrder',512, ...
            'StopbandFrequency1',.1,'PassbandFrequency1',.2,...
            'StopbandFrequency2',4.1,'PassbandFrequency2',4,...
            'SampleRate',Fs);
        delsig(:,cell)=filtfilt(delp,Fall(:,cell));
    end
    if Fs>=24
        thep=designfilt('bandpassfir','FilterOrder',512, ...
            'StopbandFrequency1',3.9,'PassbandFrequency1',4,...
            'StopbandFrequency2',12.1,'PassbandFrequency2',12,...
            'SampleRate',Fs);
        thesig(:,cell)=filtfilt(thep,Fall(:,cell));
    end
    if Fs>=30
        ap=designfilt('bandpassfir','FilterOrder',128, ...
            'CutoffFrequency1',8,'CutoffFrequency2',15,...
            'SampleRate',Fs);
        asig(:,cell)=filtfilt(ap,Fall(:,cell));
    end
    if Fs>=60
        betp=designfilt('bandpassfir','FilterOrder',128, ...
            'CutoffFrequency1',16,'CutoffFrequency2',30,...
            'SampleRate',Fs);
        betsig(:,cell)=filtfilt(betp,Fall(:,cell));
    end
    if Fs>=75
        gp=designfilt('highpassfir','FilterOrder',128, ...
            'PassbandFrequency',30,'StopbandFrequency',25,...
            'SampleRate',Fs);
        gsig(:,cell)=filtfilt(gp,Fall(:,cell));
    end
    
    hsig(:,cell)=filtfilt(hp,Fall(:,cell));
end
clear angle

    if Fs>=24
[PRTA,meanPRTA,semPRTA]=calc_PRTA(abs(hilbert(thesig)),rewards,Fs,[MouseID,'thetafreq'],saveDir,env);
[PRTA,meanPRTA,semPRTA]=calc_PRTA(angle(hilbert(thesig)),rewards,Fs,[MouseID,'thetafreq phase'],saveDir,env);
close all
    end
    
    
swin=60;
ewin=120;
rewsmall2=zeros(length(Fall),1);
num_cells=size(Fall,2);
rewratio=length(rewards)/length(Fall);
for jj=1:length(Fall)
    rewsmall2(jj)=max(rewards(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
end
[~,rewlocs]=findpeaks(rewsmall2,'MinPeakDistance',round(Fs)*2);
avgrewdist=median(diff(rewlocs));
sframe=round(swin*Fs);
eframe=round(ewin*Fs);
constrain = @(sig) (sig-min(sig(sframe:eframe)))/max(sig(sframe:eframe)-min(sig(sframe:eframe)));
for cell=1:num_cells
    gmag=abs(hilbert(gsig(:,cell)));
    amag=abs(hilbert(asig(:,cell)));
    
    %         for n=2:2
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    %     offset=max(Fall(:,cell))/max(Fall(sframe:eframe,cell));
    offset=1;
    suptitle(['Oscillations for cell: ',num2str(cell), ' Mouse: ',MouseID])
    a1=subplot(3,1,1);
    hold on
    %     plot(times,Fall(:,cell)/max(Fall(sframe:eframe,cell))+offset*2,...
    %         times,hsig(:,cell)/max(hsig(sframe:eframe,cell))+offset,...
    %         times,lsig0(:,cell)/max(lsig0(sframe:eframe,cell)));
    %     if Fs==100
    %         plot(times,gsig(:,cell)/max(gsig(sframe:eframe,cell))+3*offset);
    %     end
    %     plot(times,pos/max(pos)+min(Fall(sframe:eframe,1)-1));%,times,(Fall(:,1)))
    plot(times,constrain(Fall(:,cell))+offset*2,...
        times,constrain(hsig(:,cell))+offset,...
        times,constrain(lsig(:,cell)));
    legend('dF/F','Above .05','Below .05')
    if Fs==100
        %plot(times,constrain(gsig(:,cell))+3*offset);
        plot(times,constrain(gmag)+3*offset);
        plot(times,constrain(amag)+4*offset);
        legend('dF/F','Above .05','Below .05','Gamma','Beta')
    end
    plot(times,pos(:,cell)/max(pos(:,cell))-1);%,times,(Fall(:,1)))
    xlim([swin ewin])
    
    title('Activity and Position')
    a2=subplot(3,1,2);
    hold on
    plot(times,constrain(Fall(:,cell))+offset*2,...
        times,constrain(hsig(:,cell))+offset,...
        times,constrain(lsig(:,cell)));
    legend('dF/F','Above .05','Below .05')
    
    if Fs==100
        plot(times,constrain(gmag)+3*offset);
        plot(times,constrain(amag)+4*offset);
        legend('dF/F','Above .05','Below .05','Gamma','Beta')
        
    end
    plot(times,speedv(:,cell)/max(speedv(:,cell))-1);%,times,(Fall(:,1)))
    title('Activity and Ball speed')
    xlim([swin ewin])
    
    a3=subplot(3,1,3);
    hold on
    plot(times,constrain(Fall(:,cell))+offset*2,...
        times,constrain(hsig(:,cell))+offset,...
        times,constrain(lsig(:,cell)));
    legend('dF/F','Above .05','Below .05')
    
    if Fs==100
        plot(times,constrain(gmag)+3*offset);
        plot(times,constrain(amag)+4*offset);
        legend('dF/F','Above .05','Below .05','Gamma','Beta')
    end
    plot(times,abs(speedy(:,cell)/max(speedy(:,cell)))-1);%,times,(Fall(:,1)))
    title('Activity and VR Speed')
    xlim([swin ewin])
    
%     subplot(3,2,2);
%     plot(fpwn,(10*log10(pwn(:,cell))))
%     
%     subplot(3,2,4)
%     plot(fpwn,(10*log10(pwn(:,cell))))
%     xlim([0 .05])
%     
%     
%     subplot(3,2,6)
%     plot(cross_corrs(round(size(cross_corrs,1)/2-Fs*2):round(size(cross_corrs,1)/2+Fs*2),cell))
%     title('Cross Corr between cell and ball speed')
%     linkaxes([a1,a2,a3],'x');
%     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'close look ', num2str(cell),',', env,'.jpg']);
%     savefig([saveDir,'\',MouseID,'\', MouseID, 'close look', num2str(cell),',',env,'.fig']);
%     
%     %Plot just df/F with speed and pos and rew
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    plot(times,(constrain(Fall(:,cell))+1))
    hold on
    plot(times,speedv(:,cell)/max(speedv(:,cell)))
    plot(times,abs(speedy(:,cell))/max(speedy(:,cell))-1)
    plot(times,pos(:,cell)/max(pos(:,cell))-1)
%     plot(times,rotationvel/max(rotationvel))
    plot(times,rewsmall2+1)
    xlim([swin ewin])
    legend('dF','Ball Forward','VR Speed', 'VR Pos','Rewards')
    title(['Speed and dF for Cell ', num2str(cell)])
    
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Pos and Speed look ', num2str(cell),',',env, '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'Pos and Speed look', num2str(cell),',',env,'.fig']);
    
    
    %Plot dF VR pos and VR speed on one axis
    
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    plot(times,((Fall(:,cell))))
    hold on
    %     plot(times,speedv/max(speedv))
    plot(times,abs(speedy(:,cell))/max(speedy(:,cell)))
    plot(times,pos(:,cell)/max(pos(:,cell)))
    %     plot(times,rotationvel/max(rotationvel))
    %     plot(times,rewsmall+1)
    if ~isnan(sum(Fall(:,cell)))
    xlim([swin ewin])
    ylim([min(Fall(:,cell)), max(Fall(:,cell))])
    end
    legend('dF','VR Speed', 'VR Pos')
    title(['Speed and dF for Cell ', num2str(cell)])
    
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'VR speed same axis look ', num2str(cell),',',env,'.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'VR speed same axis look', num2str(cell),',',env,'.fig']);
    
    %plot all freq bands (env for >theta)
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    hold on
    plot(times,(constrain(Fall(:,cell))))
    plot(times,(constrain(lsig(:,cell))+1))
    plot(times,10*(rewsmall2)-1)
    plot(times,pos(:,cell)/max(pos(:,cell))-1)
    legend('dF','low freq', 'VR Pos')
        if ~isempty(Fall(:,cell))
    ylim([-1 2])
    if Fs>6
        plot(times,(constrain(delsig(:,cell))+2))
        ylim([-1 3])
        legend('dF','low freq', 'VR Pos','delta')
    end
    if Fs>=24
        plot(times,(constrain(thesig(:,cell))+3))
        ylim([-1 4])
        legend('dF','low freq', 'VR Pos','delta','theta')
    end
    if Fs>=30
        plot(times,(constrain(abs(hilbert(asig(:,cell))))+4))
        ylim([-1 5])
        legend('dF','low freq', 'VR Pos','delta','theta','alpha')
    end
    if Fs>=60
        plot(times,(constrain(abs(hilbert(betsig(:,cell))))+5))
        ylim([-1 6])
        legend('dF','low freq', 'VR Pos','delta','theta','alpha','beta')
    end
    if Fs>=75
        plot(times,(constrain(abs(hilbert(gsig(:,cell))))+6))
        ylim([-1 7])
        legend('dF','low freq', 'VR Pos','delta','theta','alpha','beta','gamma')
    end
    hold on
    %     plot(times,rotationvel/max(rotationvel))
    %     plot(times,rewsmall+1)
    xlim([swin ewin])
    title(['Freq Bands for Cell ', num2str(cell)])
        end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Freq bands, cell ', num2str(cell),',',env,'.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'Freq Bands, cell', num2str(cell),',',env,'.fig']);
    
    %         end
end    
    
    
