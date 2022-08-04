function saveAsTiff_batch_MA(nvids)
%     if nargin==1
for n=1:nvids
    [FileIn{n},PathIn{n}] = uigetfile('*.mat','Video To Convert');
    [FileName{n},PathName{n}] = uiputfile({'*.tif'},'Save Image Stack');
    FileName{n}
%     else 
%         FileName=varagin{2};
%         PathName=varargin{1};
%     end
end
    javaaddpath 'C:\Program Files\MATLAB\R2016a\java\mij.jar' 
    javaaddpath 'C:\Program Files\MATLAB\R2016a\java\ij-1.51h.jar'    
    MIJ.start;
for n=1:nvids
    load([PathIn{n},FileIn{n}]);
    if exist('chone','var')
        vid=chone;
        clear chone
    else vid=video;
        clear video
    end
                lims(1)=min(vid(:));
            lims(2)=max(vid(:));
    cmap=colormap(gray(512));
%     nParts=ceil(length(vid)/1500);
%         if FileName{n} == 0
%         else 
%             lims(1)=min(vid(:));
%             lims(2)=max(vid(:));
%             for p=1:nParts
%                 if p<nParts
%                 outLength=1500;
%                 else
%                     outLength=mod(length(vid),1500);
%                 end
%             for i=1:outLength
%                 imwrite(gray2ind(mat2gray(vid(:,:,(p-1)*1500+i),double(lims)),512),[PathName{n} FileName{n}(1:end-4) '_' num2str(p) '.tif'],...
%                     'WriteMode', 'append',  'Compression','none');
%             end
%             end
%         end
       imageJ_savefilename=strrep([PathName{n},'\',FileName{n}(1:end-4),'.tif'],'\','\\'); %ImageJ needs double slash
       imageJ_savefilename=['path=[' imageJ_savefilename ']'];
%      MIJ.createImage('chone_image', gray2ind(mat2gray(vid(:,:,((p-1)*lengthVid+1):min((p)*lengthVid,end)),double(lims)),512), true);
          MIJ.createImage('chone_image', gray2ind(mat2gray(vid,double(lims)),512), true);

        MIJ.run('Save', imageJ_savefilename);   %saves with defined filename
        MIJ.run('Close All');
end