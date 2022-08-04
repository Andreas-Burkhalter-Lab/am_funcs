function a = squeeze_angle(theta);

for i=1:length(theta)
    while (theta(i)>360) 
        theta(i) = theta(i) - 360;
    end
    while (theta(i)<0)
        theta(i) = theta(i) + 360;
    end
end
a=theta;