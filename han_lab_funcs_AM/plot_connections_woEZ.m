function [distall,corrall]=plot_connections_woEZ(Fall,ybinned,masks,cpp,varargin)
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
    corrRow=varargin{4};
    scatterRow=varargin{5};
    barRow=varargin{6};
end
if nargin>10
    age=varagin{7};
else age='none';
end
num_cells=size(Fall,2);
binsize=18;
path=smooth(ybinned);
path=(path-min(path));
pos=(path/max(path)*binsize+eps);
useinds=pos>1&pos<17;
% pos=pos(useinds);
Fall=Fall(useinds,:);

figure
title(['Correlation Matrix for mouse: ',MouseID,' no Ez'])
cors=corrcoef(Fall);
imagesc(cors)

colorbar
% set(gca,'Xtick',1:length(sorder),'XTickLabel',sorder)
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'CorrMat',env,' no Ez','.jpg']);
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
title(['Correlation with distance',' no Ez']);
xlabel('Cells')
legend('Correlation','P-Value')
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Bar Corr with Distances',' no Ez',env, '.jpg']);
figure
scatter(distall,corrall)
xlabel('Distance')
ylabel('Correlation')
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Scatter Corr with Distances',' no Ez', env, '.jpg']);
if exist('corrRow','var')
    import mlreportgen.dom.*;
    append(corrRow,TableEntry(Image([saveDir,'\',MouseID,'\',  MouseID, 'CorrMat',env,' no Ez','.jpg'])));
    append(scatterRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'Scatter Corr with Distances',' no Ez',env, '.jpg'])));
    append(barRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'Bar Corr with Distances',' no Ez', env,'.jpg'])));
    
    
end
