%%%%%%% manually assign number labels to points in an image with visibly
%%%%%%% labeled numbers

function namepoints(pointsimagein,labeled_image)

startind = 166;

pointinds = find(pointsimagein);
[pointsubsyx(:,1) pointsubsyx(:,2)] = find(pointsimagein);
npoints = length(pointinds);
pointlabels = table(NaN(npoints,1),pointsubsyx, 'VariableNames', {'name','locyx'});

if exist('pointlabels.mat','file')
    load('pointlabels.mat')
end


for ind = startind:npoints
    figstyle = get(0,'DefaultFigureWindowStyle');
    set(0,'DefaultFigureWindowStyle','docked');
    fh =     figure;
    % figure('units','normalized','outerposition',[0 0 1 1])
    image(labeled_image)
    hold all
    scatter(pointlabels.locyx(ind,2),pointlabels.locyx(ind,1),'r')
    hold off
    set(0,'DefaultFigureWindowStyle',figstyle);
    xlim(pointlabels.locyx(ind,2) + [-100 100])
    ylim(pointlabels.locyx(ind,1) + [-5 25])

    commandwindow
    inputname = input(['Point name? (ind=' num2str(ind) ') ']);
    pointlabels.name(ind) = inputname;
%     save pointlabels pointlabels
    close(fh)
end