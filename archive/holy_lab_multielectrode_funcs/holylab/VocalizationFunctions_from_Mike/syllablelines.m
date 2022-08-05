function [lines,f,t] = syllablelines(snip,threshold,nfft,samprate,scale)
%% Take some stats about the syllable for later

tcols = size(snip,2);
fcols = size(snip,1);
f = linspace(0,samprate/2,fcols);
t = (0:1:tcols - 1).*(nfft*0.5/samprate);

%% Calc the unnormalized PSDs. Maybe squaring is a good idea for findpeaks().
mag = full(abs(snip));
mag(mag==0) = threshold;
pow = mag.^2;


%% Determine peaks in each tcol. 
b = cell(1,tcols);
% Find peaks in each time bin
for i = 1:tcols
        [~,b{i}] = findpeaks(pow(:,i),'minpeakheight',scale^2,'minpeakdistance',20);
        em(i) = isempty(b{i});
end
k = find(em);
m = find(em == 0);
if ~isempty(m) % b must have some data
% Get rid of zeros
for i = 1:numel(k)
    d = abs(m-k(i));
    [~,md] = min(d);
    b{k(i)} = b{m(md(1))};
end
for i = 1:tcols
    s(i) = size(b{i},1);
end


% Calculate lines.
% Step 1. Pretend all the bs are unique lines, so make the maximum size
% matrix:
lines = zeros(sum(s),tcols);

for i = 1:tcols
    x = b{i};
    if i >= 2
        y = find(lines(:,i-1));
    end
    toremove = zeros(1,numel(x));
    for j = 1:numel(x)
        if i == 1
            lines(j,i) = x(j);
            toremove(j) = 1;
        elseif i >= 2
            d = abs(x(j) - lines(y,i-1));
            [~,k] = min(d);
            if d(k) <= 20
                lines(y(k),i) = x(j);
                toremove(j) = 1;
            end
        end
    end
    k = find(toremove);
    x(k) = [];
    if ~isempty(x)
        for j = 1:numel(x)
            lines(y(end)+j,i) = x(j);
        end
    end
end
else
    lines = [];
end

 
 %% Clean up any lines with length less with less that 2 values
 if ~isempty(lines)
   for i = 1:size(lines,1)
     k(i) = numel(find(lines(i,:)));
   end
   if ~isempty(find(k<5)) && ~isempty(find(k>=5))
   lines(find(k<5),:) = [];
 % If there are none, this is basically a syllable that's short, with one
 % line, with a bunch of power loss. Replace with simple peak frequency
 % over time.
   elseif (~isempty(find(k<5)) && isempty(find(k>=5)))
     lines = zeros(1,tcols);
     for i = 1:tcols
         if max(pow(:,i)) < scale^2
             lines(i) = 0;
         elseif max(pow(:,i)) > scale^2
             lines(i) = find(pow(:,i) == max(pow(:,i)),1);
         end
     end
   end
 elseif isempty(lines)
     lines = zeros(1,tcols);
     for i = 1:tcols
         if max(pow(:,i)) < scale^2
             lines(i) = 0;
         elseif max(pow(:,i)) > scale^2
             lines(i) = find(pow(:,i) == max(pow(:,i)),1);
         end
     end
 end
 
  z = find(lines);
 lines(z) = f(lines(z));
 
 
 %% Put lines in a structure
 if ~isempty(find(lines))
 k = size(lines,1);
 clear y;
 y(k) = struct('data',[],'cols',[],'mf',[],'slope',[],'overlap',zeros(1,numel(k)));
 
 for i = 1:k
     w = find(lines(i,1:tcols));
     y(i).data = lines(i,w(1):w(end));
     y(i).cols = w;
     y(i).mf = mean(y(i).data);
     y(i).slope = polyfit(t(w).*1000,y(i).data(y(i).data > 0),1);
     y(i).slope = y(i).slope(1);  
 end
 
 if k > 1
  for i = 1:k
     for j = 1:k
         y(i).overlap(j) = numel(intersect(y(i).cols,y(j).cols));
     end
  end
 else
     y = rmfield(y,'overlap');
 end
 clear lines;
 lines = y;
 
 else
     clear lines;
     lines = struct('data',[],'cols',[],'mf',[],'slope',[],'overlap',zeros(1,numel(k)));
 
 end
end
             
             
         