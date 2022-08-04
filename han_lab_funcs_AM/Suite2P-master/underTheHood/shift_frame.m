function out = shift_frame(I, y, x)

[Ly, Lx] = size(I);

yrange = mod([1:Ly] - y - 1, Ly) + 1;
xrange = mod([1:Lx] - x - 1, Lx) + 1;


out = I(yrange, xrange);