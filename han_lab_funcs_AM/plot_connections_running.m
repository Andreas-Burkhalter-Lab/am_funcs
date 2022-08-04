function [distall,corrall]=plot_connections_running(Fall,ybinned,Fs,masks,cpp,varargin)
if nargin>5
    saveDir=varargin{1};
    MouseID=varargin{2};
    if nargin>7
        env=varargin{3};
    else env='';
    end
else
    saveDir='F:\MA Data\Interneurons\';
    MouseID='This Mouse';
    
end
if nargin>8
    corrRow=varargin{4};
    scatterRow=varargin{5};
    barRow=varargin{6};
end
if nargin>11
    age=varargin{7};
else age='none';
end
num_cells=size(Fall,2);

dirs={'up','down'};
%% Split paths
%pos_run and neg_run are indices of respective portions of runs
[posL,run.up,negL,run.down]=split_paths(ybinned,Fs);

Fdir.up=Fall(run.up,:);
Fdir.down=Fall(run.down,:);

for ii=1:2
    figure
    title(['Correlation Matrix for mouse: ',MouseID,dirs{ii}])
    cors=corrcoef(Fdir.(dirs{ii}));
    imagesc(cors)
    
    colorbar
    % set(gca,'Xtick',1:length(sorder),'XTickLabel',sorder)
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'CorrMat',env,dirs{ii},'.jpg']);
    %compare with distance
    dists{length(cpp)}=0;
    rs{length(cpp)}=0;
    ps{length(cpp)}=0;
    rall=[];
    pall=[];
    distall=[];
    corrall=[];
    for p=1:length(cpp)
        if cpp(p)>2
            centroids=find_centroids(masks{p});
            if isequal(age,'old')
                dists{p}=pdist(bsxfun(@times,centroids,[.76,1]));
            elseif isequal(age,'new')
                dists{p}=pdist(bsxfun(@times,centroids,[.9,1.07]));
            else
                dists{p}=pdist(centroids);
            end
            %by cell
            distsquare=squareform(dists{p});
            rs{p}=zeros(cpp(p),1);
            ps{p}=rs{p};
            if p==1
                cells_in_this_plane=1:cpp(p);
            else cells_in_this_plane=(sum(cpp(1:(p-1)))+1):sum(cpp(1:p));
            end
            for cell=1:length(cells_in_this_plane)
                co=cors(cells_in_this_plane(cell),setxor(cells_in_this_plane,cells_in_this_plane(cell)));
                di=distsquare(cell,setxor(1:length(cells_in_this_plane),cell));
                corrall=[corrall,co];
                distall=[distall,di];
                if cpp(p)>3
                    [rt,pt]=corrcoef(co,di);
                    rs{p}(cell)=rt(1,2);
                    ps{p}(cell)=pt(1,2);
                end
            end
            if cpp(p)>3
                rall=[rall;rs{p}];
                pall=[pall;ps{p}];
            end
        end
    end
    %All distance
    [rsum,psum]=corrcoef(corrall,distall);
    bar([[rall;rsum(2,1)],[pall;psum(2,1)]]);
    xlab=strsplit(num2str(1:num_cells));
    xlab{num_cells+1}='All';
    set(gca,'Xtick',1:(num_cells+1),'XTickLabel',xlab)
    title(['Correlation with distance',dirs{ii}]);
    xlabel('Cells')
    legend('Correlation','P-Value')
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Bar Corr with Distances',dirs{ii}, '.jpg']);
    figure
    scatter(distall,corrall)
    xlabel('Distance')
    ylabel('Correlation')
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Scatter Corr with Distances',dirs{ii}, '.jpg']);
    
end
if exist('corrRow','var')
    import mlreportgen.dom.*;
    for ii=1:2
        append(corrRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'CorrMat',env,dirs{ii},'.jpg'])));
        append(scatterRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'Scatter Corr with Distances',dirs{ii}, '.jpg'])));
        append(barRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'Bar Corr with Distances',dirs{ii}, '.jpg'])));
        
    end
    
end