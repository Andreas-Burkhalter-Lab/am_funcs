   % plot the distances between templates
   templateFile='cycle2_berry.fine_cluster';
   templateVars=load(templateFile, '-mat');
   templates=templateVars.templates;
   channels=templateVars.channels;
   thresh=templateVars.thresh;
   medv=templateVars.medv;
   
   
   nTemplates=length(templates);
   dist=nan(nTemplates, nTemplates);
   for row=1:nTemplates
      for col=1:nTemplates
         t1=templates{row};
         t2=templates{col};
         dist(row, col)=sum((t1-t2).^2)/norm(t1)/norm(t2);
      end % for, each col
   end % for, each row
   
   figure;
   imagesc(dist, [0 4]);
   set(gca, 'ydir', 'normal');
   colorbar;
   