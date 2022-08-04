function [distall,corrall]=plot_connections_bt_planes(Fall,masks,cpp,varargin)
if nargin>3
    saveDir=varargin{1};
    MouseID=varargin{2};
    if nargin>5
        env=varargin{3};
    else env='';
    end
else
    saveDir='F:\MA Data\Interneurons\';
    MouseID='This Mouse';
    
end
if nargin>6
    part=varargin{4};
end
if nargin>7
    age=varargin{5};
else age='null';
end
if nargin>8
    etl=varargin{6};
end
num_cells=size(Fall,2);
figure
title(['Correlation Matrix for mouse: ',MouseID])
cors=corrcoef(Fall);
imagesc(cors)

colorbar
% set(gca,'Xtick',1:length(sorder),'XTickLabel',sorder)
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'CorrMat',env,'.jpg']);
%compare with distance
dists{length(cpp)}=0;
rs{length(cpp)}=0;
ps{length(cpp)}=0;
rall=[];
pall=[];
distall=[];
corrall=[];
%TODO: Fix to handle no etl case
if ~exist('etl','var');
    centroids=[];
    for p=1:length(cpp)
        if cpp(p)>2
            centroids=[centroids; find_centroids(masks{p})];
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
        else centroids=[centroids];
        end
        
    end
    centroids=[];
else
    centroids=[];
    etlref=[0 352 704 1056 1408 1760 1960 2112 2300];
    depthref=[0 14.45 44.14 73.44 110.94 155.08 164.53 177.11 181];
    [fittedetl,~]=createEtlFit(etlref,depthref);
    etliter=etl/length(cpp);
    etlvals=0:etliter:(etliter*(length(cpp)-1));
    
    if isequal(age,'old')
        etlheight=@(fittedetl) (160/181)*fittedetl(3000/2300*etl);
        for p=1:length(cpp)
            centroidstemp=find_centroids(masks{p});
            centroidstemp=[bsxfun(@times,centroidstemp,[.76,1]), ones(cpp(p),1)*max(0,etlheight(etlvals(p)))];
            centroids=[centroids;centroidstemp];
        end
    else etlheight=fittedetl;
        for p=1:length(cpp)
            centroidstemp=find_centroids(masks{p});
            if ~isempty(centroidstemp)
                if ~isnan(sum(centroidstemp))
                    centroidstemp=[bsxfun(@times,centroidstemp,[.9,1.07]), ones(cpp(p),1)*max(0,etlheight(min(etlvals(p),2300)))];
                    centroids=[centroids;centroidstemp];
                else centroidstemp=[];
                end
            else centroidstemp=[];
            end
        end
    end
end

dists=pdist(centroids);
%by cell
distsquare=squareform(dists);
rs=zeros(size(centroids,1),1);
ps=rs;
for cell=1:size(centroids,1)
    co=cors(cell,setxor(1:size(centroids,1),(cell)));
    di=distsquare(cell,setxor(1:size(centroids,1),(cell)));
    corrall=[corrall,co];
    distall=[distall,di];
    [rt,pt]=corrcoef(co,di);
    rs(cell)=rt(1,2);
    ps(cell)=pt(1,2);
end
rall=[rall;rs];
pall=[pall;ps];
% end


%All distance
if exist('etl','var')
    [rsum,psum]=corrcoef(corrall,distall);
    bar([[rall;rsum(2,1)],[pall;psum(2,1)]]);
    xlab=strsplit(num2str(1:num_cells));
    xlab{num_cells+1}='All';
    set(gca,'Xtick',1:(num_cells+1),'XTickLabel',xlab)
    title('Multiplane Correlation with distance');
    xlabel('Cells')
    legend('Correlation','P-Value')
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Multiplane Bar Corr with Distances', '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'Multiplane  Bar Corr with Distances', '.fig']);
    
    figure
    scatter(distall,corrall)
    xlabel('Distance')
    ylabel('Correlation')
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Multiplane Scatter Corr with Distances', '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'Multiplane Scatter Corr with Distances', '.fig']);
    title('Multiplane correlation with distance')
end
%By Planes
if length(cpp)>1
    for ii=1:length(cpp)
        coplanes{ii}{length(cpp)}=[];
    end
    for p1=1:length(cpp)
        if p1==1
            cells_in_p1=1:cpp(p1);
        else cells_in_p1=(sum(cpp(1:(p1-1)))+1):sum(cpp(1:p1));
        end
        figure
        hold on
        for cell=1:length(cells_in_p1)
            curcell=cells_in_p1(cell);
            for p2=1:length(cpp)
                if p2==1
                    cells_in_p2=1:cpp(p2);
                else cells_in_p2=(sum(cpp(1:(p2-1)))+1):sum(cpp(1:p2));
                end
                if ~isempty(setxor(cells_in_p2,curcell))
                    if p1==p2
                        coplanes{p1}{p2}=[coplanes{p1}{p2},cors(curcell,setxor(curcell,cells_in_p2))];
                    else
                        coplanes{p1}{p2}=[coplanes{p1}{p2},cors(curcell,(cells_in_p2))];
                    end
                end
            end
        end
        for p2=1:length(cpp)
            
            if ~isempty(coplanes{p1}{p2})
                scatter(repmat(p2,size(coplanes{p1}{p2})),coplanes{p1}{p2});
            end
            
        end
        title(['Between Planes Correlations for Plane ', num2str(p1)])
        xlabel('Plane Number')
        saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Corr by Plane',num2str(p1),env, '.jpg']);
        savefig([saveDir,'\',MouseID,'\', MouseID, 'Corr by Plane',num2str(p1),env, '.fig']);
        
        
    end
    if exist('part','var')
        import mlreportgen.dom.*;
        plTable=Table;
        plRow=TableRow;
        for ii=1:length(cpp)
            append(plRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'Corr by Plane',num2str(ii),env, '.jpg'])));
        end
        append(plTable,plRow);
        append(part,plTable);
    end
end
end