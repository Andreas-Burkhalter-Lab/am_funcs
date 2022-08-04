function T = findTargetIterative(mov)

nFramesSub = 1000;

[Ly, Lx, NT] = size(mov);
movSUB = mov(:, :, round(linspace(1, NT, nFramesSub)));
movSUB = single(reshape(movSUB, [], nFramesSub));

CC = corrcoef(movSUB);

[CCsort, isort] = sort(CC, 1, 'descend');

nFramesAvg = 20;
sumCC = sum(CCsort(1:nFramesAvg, :), 1);
[~, imax] = max(sumCC);

T = mean(movSUB(:, isort(1:nFramesAvg, imax)), 2);

T = reshape(T, Ly, Lx);

movSUB = reshape(movSUB, Ly, Lx, []);

for i = 1:10
   [yoffset, xoffset, ccmax] = crossCorrelation(T, movSUB);

   yoffset = yoffset - round(mean(yoffset));
   xoffset = xoffset - round(mean(xoffset));
   
   for j = 1:size(movSUB,3)
       % shift frames *back* by the determined amounts
       movSUB(:,:,j) = shift_frame(movSUB(:,:,j), -yoffset(j), -xoffset(j));
   end
   [~, isort] = sort(ccmax, 'descend');
   
   T = mean(movSUB(:, :, isort(1:300)), 3);
end

T = mean(movSUB, 3);

end