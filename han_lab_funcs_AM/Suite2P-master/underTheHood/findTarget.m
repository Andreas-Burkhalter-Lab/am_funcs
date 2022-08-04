function T = findTarget(mov)

[Ly, Lx, NT] = size(mov);
movSUB = mov(:, :, round(linspace(1, NT, 1000)));
movSUB = single(reshape(movSUB, [], 1000));

CC = corrcoef(movSUB);

[CCsort, isort] = sort(CC, 1, 'descend');

nFramesAvg = 10;
sumCC = sum(CCsort(1:nFramesAvg, :), 1);
[~, imax] = max(sumCC);

T = mean(movSUB(:, isort(1:nFramesAvg, imax)), 2);

T = reshape(T, Ly, Lx);
end