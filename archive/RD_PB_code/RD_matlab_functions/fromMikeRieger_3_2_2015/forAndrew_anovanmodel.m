%%
% clear;close all;clc;
load pretendData.mat;

%% 

% Ok a primer for your data:
% In this scenario, I have data for 10 neurons.
% Some are PARV some are PYR, these are a cell array of strings called
% group
% They have data for 4 "distances" I've coded as 0 0.5 1 and 2. Remember
% you're treating this as a "category" not a "measurement"

% Finally, for each member of the group, at each level of distance, I have
% a code for the neurons themselves, their IDs called neuronID.


% One useful thing is make sure your data are coded as factors:

group = nominal(group);
dist = nominal(dist);
neuronID = nominal(neuronID);


% This makes it easier to operator than as a cell array of strings for
% logical tests:

NParv = sum(group=='PARV');
NPyr = sum(group=='PYR');

% You can see I have 16 measurements associated with group id PARV and 24
% associated with PYR. Since each are 4 per distances, this is 16/4 = 4 and
% 24/4 = 6 neurons of each class. Because of the math of regression, it is
% unnecessary to have the same number of samples in each group.


% Anything coded as a nominal array takes logical tests for the String
% version of its names. dist=='0.5' and so on.

%% Model set up


% The mixed linear model solves the equation:

% y = X*b + Z*u
% Where X is called the fixed factor design matrix. These include numeric
% variables for the levels of the group and distances factors, as well as
% group X distance interaction terms.

% Z is the random effects matrix. Since you have measured each neuron
% multiple times, you get to account for random variance due to each neuron
% being idiosyncratic. Only certain neurons are in certain groups, so this
% is coded in a group-specific way. 

% A split plot or mixed design means the following:

% We will estimate the contributions of group membership, distance
% measurement, and group membership X distance measurement variances to the
% total variance observed.

% Group members are measured multiple times, so the model takes this into
% account in constructing its significance test.

% anovan() will deal directly with how to code these matrices, so you don't
% have to and honestly, it's a headache.

% The total variance sum of squares term is just sum(each data point -
% mean(all data points)^2. 

% You compute the coefficients b and u. Then you calculate yPredicted and
% determine the residual error sum((Ypredicted - Y).^2)

% To understand the effect of group membership on the linear model, you
% drop the corresponding column(s) from the design matrix X, recompute the
% coefficients, recompute the residual error. The difference between the
% error terms from the full and reduced model is the contribution of that
% factor. That is how the F test is designed.

% F = VARIANCE DUE TO FACTOR / ERROR VARIANCE
% 
% Finally, in the mixed design, the random terms in "Z" act as the
% denominator for the respective factors in which they are nested
% hierarchically. This is mathematically equivalent to collapsing distance,
% reducing the data to here N=4,6, and compute a one way test. Here
% everything is computed all at once:

% F(BETWEEN SUBJECTS, i.e. GROUP MEMBERSHIP) = VARIANCE DUE TO GROUP
% MEMBERSHIP/ VARIANCE DUE TO SUBJECT REPEATED MEASUREMENTS (Z)

% F(WITHIN SUBJECTS, i.e. DISTANCE) = VARIANCE DUE TO DISTANCE/VARIANCE IN
% REMAINING ERROR 

% F(WITHIN SUBJECTS, i.e. Group x Distance) = VARIANCE DUE TO DISTANCE
% GROUPED BY GROUP/ REMAINING ERROR

% Now I looked at your group. It looks like at many distances there is a
% value of zero. This might throw a wrench in the model, but we'll see what
% happens. You might consider a multivariate approach instead.


% Specify Grouping Variable Nesting Hierarchy. You'll notice only certain
% Neurons are in certain groups. You have to specify this design:

GroupingVars = {group,dist,neuronID};
NestingStructure = [ 0 0 0; 0 0 0; 1 0 0];

% This matrix is:

%            GROUP       DIST       NEURONID
% GROUP        0          0           0
% DIST         0          0           0
% NEURONID     1          0           0

% where the ith row is hierarchically nested inside the jth column. the
% order is specified by order of appearance in GroupingVars above

ModelStructure = [1 0 0; 0 0 1; 1 1 0; 0 1 0];

% This matrix is model terms specified down the rows, and grouping vars
% along the columns.

%                                GROUP     DIST        NEURONID

% Main Effect of Group               1       0           0
% NeuronID variance term             0       0           1
% Interaction of Group X Dist        0       1           0
% Main Effect of Dist                0       1           0


% I ordered it this way to put between subjects terms At the top and
% repeated measures terms at the bottom. NeuronID variance acts as the
% denominator for the F test in Group Membership, and the remaining error
% terms acts as the denominator for the other cases.

% An F test for NeuronID as a factor itself WILL be computed, but is not
% particularly interesting. It simply says some neurons are sig different
% from others.

% the "Split Plot" name comes from agriculture. Suppose I planted several
% crops, some to receive Treatment A and some to receive Treatment B, and
% I measured yields in those crops over time. This test leverages the fact
% that I have made multiple measurements in order to be more senstive to
% group differences by accounting for this level of crop to crop
% idiosyncracy.


% Finally, the F tests the hypothesis that at least one of the levels of
% any of these factors is different from all the others. In your case, we
% are only interested in the Main Effect of Group, and we are using this
% test to summarize its effect.

% Note the more I think about this, the more I think a multivariate
% approach is more appropriate.... See the other code.

[~,anovaTable] = anovan(y,GroupingVars,'random',3,...
                        'nested',NestingStructure,'model',...
                        ModelStructure,'varnames',{'Group','Distance','NeuronID'},'display','off');

% The 'random' flag specifies that the NeuronID term will be used to
% construct the subjects error term of the model. 'nested' specifies the
% hierarchy, and 'model' specifies the actual terms.
% I turn the display off:
% Now we will make this a little prettier:


anovaTable = anovaTable(:,[1,2,3,5,11,10,6,7]);

% This displays the output with the numerators and denominators of the F
% statistic specific. The 2nd line is the Main Effect of Group which is the
% summary statistic for the effect this has on the model.



