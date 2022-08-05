%%
colors = {'c','r','m','k','b','c','r','m','k','b','c','r','m','k','b'};
options.sliderwindow = 0;
ksnips = randsample(numel(snips),20);
for i = 1:numel(ksnips)
    snip = snips{ksnips(i)};
    [s,~,~,~,~,pfad] = syllableObjects(snip,f,thresh);
    spsngplot(snip,f,options);
    hold on;
    box off;
    for j = 1:numel(s)
        plot(find(s{j}),s{j}(find(s{j})),colors{j},'LineWidth',4);
    end
    plot(pfad,'xr');
    pause;
    close all hidden;
end
        