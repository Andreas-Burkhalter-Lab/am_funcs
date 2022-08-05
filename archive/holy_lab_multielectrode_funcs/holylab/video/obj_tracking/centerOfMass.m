function com=centerOfMass(mask)
% calc the center of mass of a mask

[rows, cols]=find(mask);
com=[mean(rows) mean(cols)];
