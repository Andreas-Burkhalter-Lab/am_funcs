function [x,iter] = snnls(A,b)
%SNNLS	Solution of sparse nonnegative least squares problems.
%       @(#)snnls.m Version 1.1 1/23/96
%       Pontus Matstoms, Linkopings Universitet.
%       Mikael Adlers, Linkopings Universitet
%       e-mail: pomat@math.liu.se, maiadl@mai.liu.se
%
%       x=snnls(A,b) solves the constrained linear least squares problem
%
%       min_x ||Ax-b||_2  s.t. x >= 0
%
%       Numerical method: Block principal pivoting (Portugal et al. (1993))
%
%       Matlab version 4 and sqr (Matstoms (1994)) assumed.

% Check the input data.

if nargin ~= 2,
  disp('??? Error using ==> snnls')
  disp('Too many input arguments.')
end

% Minimum degree analysis of A.

%Pc=colmmd(A);
Pc = colamd(A);
A=A(:,Pc);

% Initiate some variables.

[m,n]=size(A);

F=[];  G=1:n;                     % x(F) and Y(G) basic variables.
x=0;   y=-A'*b;
p=3;
nf=n+1;
iter=0;

% Start of main loop

while any([x(F)<0; y(G)<0]),
  iter=iter+1;
  H1=F(find(x(F)<0));
  H2=G(find(y(G)<0));
  H12=sort([H1; H2]);             % H12 set of infeasibilities
  nf0=length(H12);
  if nf0 < nf,                    % Does |H12| decrease?
    nf=nf0;
    p=3;
  else
    if p>=1,                      % Try atmost p times to get decreasing nf.
      p=p-1;
    else
      r=H12(nf0);                 % Use Murty's method until decrease.
      if find(H1 == r),
        H1=r;
        H2=[];
      else
        H1=[];
        H2=r;
      end
    end
  end
  
% Update F and G.

  FLAG=zeros(n,1); FLAG([F; H2])=ones(size([F; H2])); FLAG(H1)=zeros(size(H1));
  F=find(FLAG);
  
  FLAG=zeros(n,1); FLAG([G; H1])=ones(size([G; H1])); FLAG(H2)=zeros(size(H2));
  G=find(FLAG);
  
% Compute x(F) and y(G).

  %x=zeros(n,1); R=chol(A(:,F)'*A(:,F)); x(F)=R\(R'\(A(:,F)'*b));
  x=zeros(n,1);
  [R,ptmp]=chol(A(:,F)'*A(:,F));
  x(F)=R\(R'\(A(:,F)'*b));
  y=zeros(n,1); y(G)=A(:,G)'*(A(:,F)*x(F)-b);
  
end   % end while


x(F)=A(:,F)\b;
x(Pc)=x;

end   % end snnls

	
      

