function shufx = shuffle(x,t)

% This function takes a vector x and randomly shuffles it t times
% to finally store it as shufx


for i=1:t
	shufx = x;
	shufx(2,1:length(x)) = rand(1,length(x));
	shufx = shufx';
	shufx = sortrows(shufx,2);
	shufx = shufx';
	shufx = shufx(1,:);
	
end



end
