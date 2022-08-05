    imold = double(imread('cameraman.tif'));
    im=repmat(imold, 4, 4*16);
    
    [X1,X2] = ndgrid(1:size(im,1),1:size(im,2));
    g = X1; g(:,:,2) = X2;
    g = g + randn(size(g));
    % Do the interpolation
    imi_st= imqinterp_st(g,im);
    imi_mt= imqinterp_mt(g,im);
    % figure; subplot(1,2,1); imshowsc(im); subplot(1,2,2); imshowsc(imi);
    if(~isequalwithequalnans(imi_st, imi_mt))
       error('mt doesn''t match non-mt\n'); 
    else
       fprintf('----passed self test-------\n');
    end
    
    fprintf('meature time of non-mt version:\n');
    tic; for i=1:5; imi_st = imqinterp_st(g,im); end; duration_st=toc
    fprintf('measure time of mt version: \n');
    tic; for i=1:5; imi_mt = imqinterp_mt(g,im); end; duration_mt=toc
    fprintf('the speedup is: %.2f\n', duration_st/duration_mt);
    