%PlotTwoHists.m: This function takes in two vectors of numbers (of generally different)
%	lengths, and plots them together on a single histogram.  This is useful, for example,
%	for plotting spike counts sorted by choice or disparity, etc.  GCD, 8/10/00
function PlotTwoHists(dist1, dist2)

%dist1 and dist2 need to be column vectors to use 'hist', so check for this
if ( (size(dist1,2)>1) & (size(dist1,1)==1) )	% a row vector, so transpose
   dist1 = dist1';
end
if ( (size(dist2,2)>1) & (size(dist2,1)==1) )	% a row vector, so transpose
   dist2 = dist2';
end

%find maximum length of vectors and use this to fill Mtemp with NaN's
max_len = max(length(dist1), length(dist2));
Mtemp = ones(max_len, 2)*NaN;

%now, put the vectors into Mtemp
Mtemp(1:length(dist1), 1) = dist1;
Mtemp(1:length(dist2), 2) = dist2;

%now, plot the histogram
hist(Mtemp);

return;