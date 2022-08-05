function winnorm = windownorm(window,noverlap)
% WINDOWNORM: compute the normalization factor for gluing windows together
% winnorm = windownorm(window,noverlap)
winlen = length(window);
shift = winlen - noverlap;
if (shift < 0 | shift > winlen)
        error('Windows must overlap properly');
end
winnorm = window;
nshift = ceil(winlen/shift);
forward = window;
backward = window;
for k = 1:nshift
        forward(1+shift:end) = forward(1:end-shift);
        forward(1:shift) = 0;
        backward(1:end-shift) = backward(1+shift:end);
        backward(end-shift+1:end) = 0;
        winnorm = winnorm + forward + backward;
end
