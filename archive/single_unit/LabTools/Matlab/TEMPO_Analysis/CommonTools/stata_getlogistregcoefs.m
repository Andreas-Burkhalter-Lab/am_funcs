function [coeffs] = Stata_getLogistRegCoefs(logfile)
% Pokes through the Stata log file from the logistic regression
% analysis ('blogit') and returns a matrix of parameter estimates and p-values.


%read in the Stata log file
[s, w] = unix(['type ' logfile]);

%find all occurences of the vertical bar
wherebars = find(w=='|');

% starting at the fourth bar, read in the data, a line at a time
% the parameters here are highly specific to the format of the log file
% Note: If you change the format of the log file, you will need to change this
len = (wherebars(5)-9) - (wherebars(4)+1); 	%length of data line to get
coeffs = [];
for i=4:size(wherebars,2)
   str = w(1, wherebars(i)+1:wherebars(i)+len);
   [temp] = sscanf(str, '%f');
   coeffs = [coeffs; temp'];
end

