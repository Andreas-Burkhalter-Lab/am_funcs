function snipplot(snips,twhis,samprate,nfreq,output)

dur = diff(twhis);
f = linspace(0,samprate/2,nfreq+1);
options.sliderwindow = 0;
for i = 32:length(snips)
    t = linspace(0,dur(i),size(snips{i},2));
    spsngplot(snips{i},f,t,options);
    title(['Assigned type: ',output.names{i}]);
    print('-depsc2',[sprintf('%03d',i),'.eps']);
    close;
end

end

