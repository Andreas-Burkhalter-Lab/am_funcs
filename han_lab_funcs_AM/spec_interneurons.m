function spec_interneurons(Fall,Fs,pos,speedv,varargin)
if nargin>4
    saveDir=varargin{1};
    MouseID=varargin{2};
    if nargin>6
        env=varargin{3};
    else env='';
    end
else
    saveDir='F:\MA Data\Interneurons\';
    MouseID='This Mouse';
    
end
if nargin>7
    part=varargin{4};
end
params.fpass = [0 Fs/2];
params.Fs=Fs;
params.tapers=[5 9];
movingwin = [5 1];
num_cells=size(Fall,2);
times=linspace(1,length(Fall)/Fs,length(Fall));

for jj=1:num_cells
    if ~isnan(sum(Fall(:,jj)))
        % jj=1;
        %Spectrogram
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        suptitle(['Spectrogram for cell: ',num2str(jj), ' Mouse: ',MouseID])
        
        a1=subplot(4,1,1:2);
        [specgram_data.s,specgram_data.t,specgram_data.f]=mtspecgramc(Fall(:,jj),movingwin,params);
        imagesc(specgram_data.t,specgram_data.f, log(specgram_data.s'));
        
        
        a2=subplot(4,1,3:4);
        hold on
        jj
        size(Fall)
        size(pos)
        plot(times,pos(1:length(Fall),1)/max(pos(1:length(Fall),1))+(max(Fall(:,jj))))
        
        plot(times,speedv(:,1)/max(speedv(:,1))+(min(Fall(:,jj))-1))
        plot(times,Fall(:,jj));
        xlim([0 max(times)])
        linkaxes([a1,a2],'x');
        saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'spectrogram Tres', num2str(jj), env,'.jpg']);
        savefig([saveDir,'\',MouseID,'\', MouseID, 'spectrogram Tres', num2str(jj),env,'.fig']);
    end
end
movingwin = [50 20];
for jj=1:num_cells
    % jj=1;
    %Spectrogram
    if ~isnan(sum(Fall(:,jj)))
        
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        suptitle(['Spectrogram for cell: ',num2str(jj), ' Mouse: ',MouseID])
        
        a1=subplot(4,1,1:2);
        [specgram_data.s,specgram_data.t,specgram_data.f]=mtspecgramc(Fall(:,jj),movingwin,params);
        imagesc(specgram_data.t,specgram_data.f, log(specgram_data.s'))
        
        
        a2=subplot(4,1,3:4);
        hold on
        plot(times,pos(:,1)/max(pos(:,1))+(max(Fall(:,jj))))
        
        plot(times,speedv(:,1)/max(speedv(:,1))+(min(Fall(:,jj))-1))
        plot(times,Fall(:,jj));
        xlim([0 max(times)])
        linkaxes([a1,a2],'x');
        saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'spectrogram Fres', num2str(jj),env, '.jpg']);
        savefig([saveDir,'\',MouseID,'\', MouseID, 'spectrogram Fres', num2str(jj),env,'.fig']);
    end
end
if exist('part','var')
    import mlreportgen.dom.*;
        specTable=Table;
        specRow=TableRow;
    for jj=1:num_cells
            if ~isnan(sum(Fall(:,jj)))

        if Fs>10 %tres
            append(specRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'spectrogram Tres', num2str(jj), env,'.jpg'])));
        else
            append(specRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'spectrogram Fres', num2str(jj), env,'.jpg'])));
        end
            end
    end
    append(specTable,specRow);
    append(part,specTable);
end
