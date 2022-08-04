%%
clear;close all;clc;
load pretendData.mat;

%% 

% Ok, in this scenario, we're going to do something a little different.
% We're actually going to compute this in a multivariate way.

% First, I'm going to code the variables as before as categorical arrays.
% I'm going to leave distance numeric so that it sorts in order and doesn't
% do anything funky with how the distances are ordered.

group = nominal(group);
neuronID = nominal(neuronID);

% Now, I'm going to summarize these results in tabular form:

Y = dataset(y,group,dist,neuronID,'VarNames',{'y','group','dist','neuronID'});

% Now I'm going to sort, first by distance, then by ID, so that my matrix
% is ready to reshape. 

Y = sortrows(Y,{'dist','neuronID'});

% Now this is easy. I have 10 IDs, so I'm going to reshape y so that
% neurons are down the rows, and by-distance measurements are across the
% columns:


yy = reshape(Y.y,length(unique(Y.neuronID)),length(unique(Y.dist)));

% You can see yy is a 10 x 4. I have 10 neurons.

% Now I need the corresponding group labels. My trick is to reshape group
% and then keep the first row, and keep just the first column.

gg = reshape(Y.group,length(unique(Y.neuronID)),length(unique(Y.dist)));

% As a quick QC you should be sure each column is identical.

gg = gg(:,1);

% Great, now we're ready to go.

%% Model set up

% Ok, so multivariate ANOVA asks the simple question, is at least one group
% different from the others. But instead of predicting a vector, we are
% asked to predict a matrix M x N Y where rows are observations, and
% columns are different dependent variables, in this case we are
% considering normalized signal at each level of distance to be its own
% type of measurement.

% Multivariate ANOVA using multidimensional scaling, to reduce the data
% into what is something like principal component analysis. It asks the
% ANOVA question along each dimension of the projected data and group
% membership may make a difference on one or more dimensions.

% if d = 0, then the means are not different from one another on any
% dimension

% if d = 1, then the means are different along a single dimension ( a line)
% if d = 2, then the means are different in a plane
% if d = 3, then the means are different in a 3D space and so on


[d,p,stats] = manova1(yy,gg);

% Ok, so multivariate stats. To explore your date, compute the projection
% of your data onto the new space:

yy2 = yy*stats.eigenvec;

% If you look at the eigenvectors themselves, the largest ranking number in
% any column corresponds to the variable which loads the heaviest on that
% column for the projection.

% stats.eigenvec(1,1) is 30.22 this means Distance #1 is most highly
% correlated with the axis described by the new column in my new reduced
% data space.

% Each column of eigenvec is the linear combination of dependent variables
% to create a new component. To understand whether or not we have reduced
% this space from our original multivariate space, look at the eigenvalues.

PDE = stats.eigenval./sum(stats.eigenval);

% This is the proportion of group discrimination explained by each variable
% in the new reduced space. The MANOVA is like MDS - it searches for a
% reduced space that maximizes the variance that best discriminates the
% Group Membership. plot(PDE) is called  Scree plot. If you have similar
% variance explained across all variables, then you haven't really done
% ANYTHING by computing the new space. 

% In my case, 100% of the variance that best discriminates my groups can be
% explained by one of my variables. That's not surprising because my
% variables are actually random numbers :)


% To explore the projections, let's do gscatter. Suppose we had D = 0, p <
% 0.05, but D = 1, p > 0.05.

% That would imply that our discrimination was observed at the level of a
% single dimension. Scatter the projection yy2 across the first column:

gscatter(yy2(:,1),zeros(length(gg),1), gg)

% I've added a column of zeros because gscatter wants 2 for plotting. Ok,
% so I've put everyone on a line, and you can see that MANOVA has done a
% good job of distinguishing the two groups, however, the within group
% variance is such that in my actual model the groups are not sig
% different.

% If I had D = 0, p < 0.05, D = 1, p < 0.05, & D = 2 p > 0.05, this implies
% that the groups are different, and probably reside in the plane made by
% the first 2 projection variables.:

gscatter(yy2(:,1),yy2(:,2),gg)

% If you look in the plane, you can see again MANOVA did a good job of
% discriminating my groups but still I know I have no sig differences
% because D = 0 and p > 0.05, meaning there is not sufficient difference
% across any dimension of my reduced space.

% If D = 0,1,2 p <0.05 but D = 3 p> 0.05 this means the group difference
% likely resides in a 3D space, and so on.


% Once you get above 3D space, you can't plot things anymore. Instead we
% turn to Euclidean distances. To understand this, let me add the MEANS of
% each groups to the scatter plot:

gscatter(yy2(:,1),yy2(:,2),gg);
mns = grpstats(yy2,gg);
hold on;
plot(mns(1,1),mns(1,2),'dr')
plot(mns(2,1),mns(2,2),'dc')

% Ok great. To understand why I have no sig differences in my data, or do,
% I can understand the distance in the plan from each of my data points to
% its group mean. If the distance between the group means in the plane made
% by the first two eigenvectors is bigger than the distance between any one
% data point to its own mean, that means I have probably a sig difference
% in my groups. 

% stat.mdist is the Euclidean distance between each datapoint across all
% dimensions and its group mean. Actually,  its the Mahalanobis distance,
% which is the Euclidean distance weighted by covariance in the data. It's
% better than Euclidean distance because it takes into account any
% dependence of one of the variables on another.

boxplot(stats.mdist,gg)

% In this case you can see these distances are not particular different
% from each other. That's sort of a test of homogeneity of variance: Each
% group is similarly distanced from its own mean. So even if my differences
% were sig I might expect this.

plot(repmat(1,sum(gg=='PARV'),1),stats.mdist(gg=='PARV'),'dk');
hold on;
plot(repmat(2,sum(gg=='PYR'),1),stats.mdist(gg=='PYR'),'dk');
plot([0 3],[stats.gmdist(1,2),stats.gmdist(1,2)]);
set(gca,'Xtick',[1,2]);
set(gca,'Xticklabel',{'PARV','PYR'});

% One thing that you will appreciate is that that two group means
% (represented by the line) are separated by a distance which is for at
% least one group well within the distribution of any of that groups
% members from its own group mean. And for the other group its pretty
% close.

% This is not surprising because I know D = 0 and p > 0.05, so my groups
% are NOT sig different, HOWEVER if they were, this would allow me to
% quantify the magnitude of difference between them in a multivariate way.


