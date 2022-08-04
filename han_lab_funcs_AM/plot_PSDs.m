function plot_PSDs(Fall,Fs,varargin)
if nargin>2
    saveDir=varargin{1};
    MouseID=varargin{2};
    if nargin>4
        env=varargin{3};
    else env='';
    end
else
    
    saveDir='F:\MA Data\Interneurons\';
    MouseID='This Mouse';
end
if nargin>5
    part=varargin{4};
    disp('Found part')
end
win=hanning(2^(floor(log2(length(Fall)))-1));
nover=2^(floor(log2(length(Fall)))-2);
nff=2^(floor(log2(length(Fall)))-1);
num_cells=size(Fall,2);
%% power spectrum
[pwn,fpwn]=pwelch(Fall(:,:),win,nover,nff,Fs);

%% plotting
num_square=9;
num_figs=ceil(num_cells/num_square);
for f=1:num_figs
    fig_start=(f-1)*num_square+1;
    if f<num_figs
        fig_end=fig_start+(num_square-1);
    else fig_end = num_cells;
    end
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    suptitle(['PSD for cells: ',num2str(fig_start),' to ',num2str(fig_end), ' Mouse: ',MouseID])
    
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        hold on
        %         len=round(length(fpwn)/10);
        plot(fpwn,(10*log10(pwn(:,cell))))
        %         if ismember(cell,npil)
        %             title('Neurites')
        %         elseif ismember(cell, nothing)
        %             title('Nothing Visible')
        %         elseif ismember(cell,missing)
        %             title(['(Missing) Cell ', num2str(cell)])
        %         else
        %             title(['Cell ', num2str(cell)])
        %         end
        title(['Cell ', num2str(cell)])
        
        splot=splot+1;
        
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'PSD', num2str(f), env,'.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'PSD', num2str(f),env,'.fig']);
    
end
if exist('part','var')
    import mlreportgen.dom.*;
    pTable=Table;
    pRow=TableRow;
    for f=1:num_figs
        append(pRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'PSD', num2str(f), env,'.jpg'])));
    end
    append(pTable,pRow);
    append(part,pTable);
    
end
