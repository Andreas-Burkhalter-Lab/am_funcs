function [lambda,u,v,err,imtrx] = prodcurrent(K,w,ctensor)
% PRODCURRENT: approximate mitral cell current as a product
%
% The mitral cell current is calculated over a matrix of concentrations
% whose paramterization is determined by an input tensor.  The current is
% then approximated as Imtrx = lambda * u*v', where u and v are column
% vectors.
%
% This functional form may be useful for designing circuits which satisfy
% certain properties, e.g., concentration invariance.
%
% Usage:
%   [lambda,u,v,err] = prodcurrent(K,w,ctensor)
% where
%   K is a ncells-by-2 matrix of association constants
%   w is a vector of length ncells giving synaptic weights
%   ctensor is a n-by-m-by-2 tensor, where c(i,j,1) is the concentration
%     of compound #1 at vertex (i,j), and c(i,j,2) is the concentration
%     of compound #2 at this vertex.  For example, if you want to
%     decompose the mitral cell current as a product of f(c1/c2) *
%     g(c1*c2), then ctensor(i,:,:) should correspond to lines of
%     constant c1/c2, and ctensor(:,j,:) should correspond to lines of
%     constant c1*c2.
%     If at any vertex, the sum of concentrations is larger than 1, this
%     vertex is not used for generating the best-fit to Imtrx. 
% and
%   lambda,u,v are described above
%   err is the mean-square error 
%     ( i.e., sum(sum((Imtrx - lambda*v*v').^2)) )
%
% [lambda,u,v,err,imtrx] = prodcurrent(...) also returns the original
%   current matrix.  Any concentration entries which fail to be sensible
%   (see explanation of ctensor above) will have 0 current.  (This is the
%   right thing to do for the SVD algorithm.)
%
% See also: CTENSRATIO.

% Tim Holy, 01-21-2002.  
[n,m,p] = size(ctensor);
if (p ~= 2)
  error('The concentration tensor must be n-by-m-by-2');
end
for i = 1:n
  for j = 1:m
    % Check to see if the concentration is valid
    if (sum(ctensor(i,j,:)) > 1)
      imtrx(i,j) = 0;
    else
      imtrx(i,j) = calccurrent(K,w,ctensor(i,j,:));
    end
  end
end  
[U,S,V] = svd(imtrx);
dS = diag(S);
lambda = dS(1);
u = U(:,1);
v = V(:,1);
err = sum(dS(2:end).^2);  % The sum of the squares of the remaining
                          % singular values is equal to the error.
