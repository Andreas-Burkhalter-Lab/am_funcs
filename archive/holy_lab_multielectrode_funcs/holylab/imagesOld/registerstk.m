function registerstk(options)
  if (nargin < 1)
    options = struct;
  end
  if ~isfield(options,'images')
    if ~isfield(options,'filein')
      [filein,pathin] = uigetfile({'*.stk';'*.tif'});
      if isempty(filein)
        return;
      end
      filein = [pathin filein];
    end
    [pathstr,basename,extname] = fileparts(filein);
    if strcmp(extname,'.stk')
        [s,nimg] = tiffread(filein);
    elseif strcmp(extname,'.tif')
      tNumToDiscard=2;
      iminfs =  imfinfo(filein);
      iminfs = iminfs(tNumToDiscard+1:end); % jason: discard the first tNumToDiscard frames 4831
      nimg = length(iminfs);
      for i = 1:nimg
          s(i).data = imread(filein,i+tNumToDiscard);
          s(i).width = iminfs(i).Width;
          s(i).height = iminfs(i).Height;
      end
    end 
  else
    s = options.images;
    nimg = length(s);
  end
  % Average some frames together
  if ~isfield(options,'avgrange')
    options.avgrange = [1 30];
  end
  imavg = zeros(size(s(1).data));
  for i = options.avgrange(1):options.avgrange(end)
    imavg = imavg + double(s(i).data);
  end
  imavg = imavg/(diff(options.avgrange([1 end]))+1);
  % Pick features in this average image
  figure
  imagesc(imavg)
  colormap(gray)
  [x,y] = ginput(6);
  ptsbase = [x,y];
  cpcorrsize = 10;
  % Now loop over images, using correlation to find the points in each
  % image corresponding to the ones in the averaged image
  ptscur = ptsbase;
  figure('Position',[4         623        1259         283]);
  subplot(1,4,1);
  imagesc(imavg);
  colormap(gray);
  for i = 1:size(ptscur,1)
    rectangle('position',[ptscur(i,1:2) 0 0] + [0 0 1 1]*2*cpcorrsize+1 - [1 1 0 0]*cpcorrsize,'EraseMode','xor');
  end
  for i = 1:nimg
    % Calculate the shift
    ptscur = mycpcorr(ptscur,ptsbase,s(i).data,imavg,cpcorrsize);
    ptscur - ptsbase
    % Now shift images
    i
    tform = cp2tform(ptscur,ptsbase,'linear conformal');
    B = imtransform(s(i).data,tform,'XData',[1 size(s(i).data,2)],'YData',[1 size(s(i).data,1)]);
    subplot(1,4,2)
    imagesc(s(i).data);
    for j = 1:size(ptscur,1)
      rectangle('position',[ptscur(j,1:2) 0 0] + [0 0 1 1]*2*cpcorrsize+1 - [1 1 0 0]*cpcorrsize,'EraseMode','xor');
    end
    subplot(1,4,3);
    imagesc(B);
    subplot(1,4,4);
    imagesc(double(B)-imavg);
    s(i).data = B;
    %pause
    drawnow
  end
  % Save new images
  if ~isfield(options,'fileout')
    [fileout,pathout] = uiputfile;
    if isequal(fileout,0)
      return;
    end
    options.fileout = [pathout fileout];
  end
  %save(options.fileout,'s');
  imwrite(s(1).data,options.fileout,'tif', ...
          'writemode','overwrite','compression','none')
  for i = 2:length(s)
    imwrite(s(i).data,options.fileout,'tif',...
            'writemode','append','compression','none');
  end