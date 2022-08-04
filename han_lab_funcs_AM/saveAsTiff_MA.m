function saveAsTiff_MA(vid,varargin)
    if nargin==1
    [FileName,PathName] = uiputfile({'*.tif'},'Save Image Stack');
    else 
        FileName=varagin{2};
        PathName=varargin{1};
    end
    
javaaddpath 'C:\Program Files\MATLAB\R2017a\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2017a\java\ij-1.51n.jar'
    cmap=colormap(gray(512));
    lengthVid=1000;
    nParts=ceil(length(vid)/lengthVid);
        if FileName == 0
        else 
            lims(1)=min(vid(:));
            lims(2)=max(vid(:));
            MIJ.start;    %calls Fiji
%             for p=1:nParts
%                 if p<nParts
%                 outLength=lengthVid;
%                 else
%                     outLength=mod(length(vid),lengthVid);
%                 end
%             for i=1:outLength
%                 imwrite(gray2ind(mat2gray(vid(:,:,(p-1)*lengthVid+i),double(lims)),512),[PathName FileName(1:end-4) '_' num2str(p) '.tif'],...
%                     'WriteMode', 'append',  'Compression','none');
%             end
%             
%                 opt.big='true';
%                 opt.overwrite='true';
%                 saveastiff(vid,[PathName FileName(1:end-4) '_' num2str(p) '.tif'],opt);
%             final_filename=[tiffpath new_tifffilename]; %files need to have path
%            imageJ_savefilename=strrep([PathName,'\',FileName(1:end-4),num2str(p),'.tif'],'\','\\'); %ImageJ needs double slash
           imageJ_savefilename=strrep([PathName,'\',FileName(1:end-4),'.tif'],'\','\\'); %ImageJ needs double slash
       imageJ_savefilename=['path=[' imageJ_savefilename ']'];
%      MIJ.createImage('chone_image', gray2ind(mat2gray(vid(:,:,((p-1)*lengthVid+1):min((p)*lengthVid,end)),double(lims)),512), true);
          MIJ.createImage('chone_image', gray2ind(mat2gray(vid,double(lims)),512), true);

        MIJ.run('Save', imageJ_savefilename);   %saves with defined filename
        MIJ.run('Close All');
%         MIJ.run('Collect Garbage');
%         pause(1);
%         MIJ.exit;
%             end
        end