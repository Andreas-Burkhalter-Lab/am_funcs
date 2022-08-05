th = 30*pi/180;
%im = ones(128,128,'uint8');
im = ones(128,128);
x = 1:size(im,2);
y = round(x * tan(th));
indx = sub2ind(size(im),y,x);
im(indx) = 2;
figure; imagesc(im)

imdec = im(:,1:2:end);
imdec = reshape(imdec,[1 size(imdec)]);
ops = struct('pixel_spacing',[1 1 2]);
angles = [th 0 0];

[imr,coords] = rotatestack(imdec,angles,ops);
figure; imagesc(squeeze(imr))
imrf = imr; imrf(isnan(imr)) = -1;
figure; imagesc(squeeze(imrf))