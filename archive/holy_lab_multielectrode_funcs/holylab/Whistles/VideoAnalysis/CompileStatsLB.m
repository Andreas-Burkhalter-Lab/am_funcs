function CompileStatsLB(vectv,behavnames,fnames,statnames)
nmice = size(vectv,1);
nfiles = size(vectv,2);
nbehavs = size(vectv,3);
nstats = size(vectv,4);
smean = zeros(nfiles,nbehavs,nstats);
serr = means;
for i = 1:nfiles
        for j = 1:nbehavs
                for k = 1:nstats
                        smean(i,j,k) = mean(vectv(:,i,j,k);
                        serr(i,j,k) = std(vectv(:,i,j,k))/sqrt(nmice);
                end
        end
end
dims = CalcSubplotDims(nfiles);
for k = 1:nstats
        figure
        for i = 1:nfiles
                subplot(dims(1),dims(2),i)
                errorbar(1:nstats,smean(i,:,k),serr(i,:,k))
                set(gca,'XTickLabel',behavnames);
                title(fnames{i});
        end
        suptitle(statnames{k})
end
