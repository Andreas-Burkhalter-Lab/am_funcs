%AngleWrap.m: this is a little function to take an angle (in degrees) and
%wrap it around into various ranges: 360, 180, and 90 degrees.
function [ang360, ang180, ang90] = AngleWrap(angle)

%first contrain the angle into the range of 0-360 deg
ang360 = angle;
while (ang360 > 360)
    ang360 = ang360 - 360;
end
while (ang360 < 0)
    ang360 = ang360 + 360;
end

%now, fold around into 180deg range
if (ang360 > 180)
    ang180 = abs(ang360 - 360);
else
    ang180 = ang360;
end

%now fold into 90 deg
ang90 = ang180 - 90;
ang90 = 90 - abs(ang90);

return;