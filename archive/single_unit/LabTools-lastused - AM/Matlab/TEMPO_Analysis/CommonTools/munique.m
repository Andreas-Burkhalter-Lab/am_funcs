function [D, c] = munique(R);
%$$$ Congratulations to Richard Aufrichtig who won this month's M-file
%contest  with the following elegant and blazingly fast solution:  The
%matrix, D, contains the unique rows of R. The vector, c, contains a count
% of the number of occurrences of each row of D in the original matrix,
%R.    The first line in Richard's solution simplifies the problem of
%finding unique  rows in a matrix to the problem of finding unique
%elements in a vector. He  does this by multiplying R by a random vector
%with length equal to the number  of columns of R.  
%  
%  Theoretically, of course, there are an infinity of rows that would multiply
%  the random vector to give the same element. But the probability that any two
%  of these rows would be in the same matrix is vanishingly small. So, the
%  solution works.

Q = R*rand(size(R,2),1);
[y, i]=sort(Q);
s=[1; diff(y)];
  
D = R(i(s~=0),:);
c=diff(find([s ;1]~=0));


 
