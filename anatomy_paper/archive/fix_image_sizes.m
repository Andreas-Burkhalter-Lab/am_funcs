% enlarge  small images
for i = 5:14
    imname = ['17185_s' num2str(i) 'b_xmicro3_2_cell-dense_roi.png'];
    imorig = loadbw(imname);
    imresized = imresize(imorig,[1040,1392]);
    delete(imname)
    savetif(imresized,getfname(imname))
end
    