%%% find the elevation and azimuth of the screen by inputing the distance
%%% of the screen center above and temporal of the eye
% arg1 = distance temporal of eye in cm
% arg2 = distance above eye in cm
%%% displays azimuth, elevation in degrees

%%% last edited AM 6/30/16

function getAngle(distUp, distTemporal)

viewDist = 16; % distance from eye to screen center in cm

az = rad2deg(asin(distUp/viewDist));
el = rad2deg(asin(distTemporal/viewDist));

fprintf('Azimuth = %g deg\n',az)
fprintf('Elevation = %g deg\n',el)
fprintf('View distance = %gcm\n',viewDist)

end