% halfrect.m -- half-wave rectifies sinusoidal data for curve fitting

function rectdata = halfrect(ydata)

rectdata = ydata;
rectdata(rectdata < 0) = 0;

return;