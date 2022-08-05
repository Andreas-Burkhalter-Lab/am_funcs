function result=color2gray(rgb)
   rgb=double(rgb);
   result=(rgb(:,:,1)*11+rgb(:,:,2)*16+rgb(:,:,3)*5)/32;
   result=result/255;
   