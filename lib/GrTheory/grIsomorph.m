function [IsIsomorph,Permut]=grIsomorph(E1,E2,n)
% Function [IsIsomorph,Permut]=grIsomorph(E1,E2,n) 
% compare two simple graphs (without loops and multiple edges) 
% about their isomorphism.
% Input parameters: 
%   E1(m1,2) - the edges of 1st graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m1 - number of edges;
%   E2(m2,2) - the edges of 2nd graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m2 - number of edges;
%   n - number of vertexes (optional); if this parameter 
%     are omitted, by default n1=max(E1(:)); n2=max(E2(:));
%     if n1~=n2, the graphs are nonisomorphic.
% Output parameters:
%   IsIsomorph = true for isomorphic graphs;
%   IsIsomorph = false otherwise;
%   For nonisomorphic graphs Permut = [];
%   For isomorphic graphs Permut(1:n,1) is permutation 
%   from vertexes numbers G1 to vertexes numbers of G2.
% Needed other products: ZerOne Toolbox.
% This software may be free downloaded from site:
% http://www.apmath.spbu.ru/grafomann/download.html
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<2,
  error('There are no input data!')
end
[m1,n1,E1] = grValidation(E1); % E1 data validation
[m2,n2,E2] = grValidation(E2); % E2 data validation
IsIsomorph=false;
Permut=[];
if nargin>=3,
  n1=n;
  n2=n;
end
if (m1~=m2)|(n1~=n2),
  return;
end
A1=zeros(n1);
A1((E1(:,2)-1)*n1+E1(:,1))=1;
A1=A1+A1'; % adjacency matrix A1
A2=zeros(n2);
A2((E2(:,2)-1)*n2+E2(:,1))=1;
A2=A2+A2'; % adjacency matrix A2
[ef,P]=psimilar(A1,A2);
IsIsomorph=(ef>0);
if IsIsomorph,
  [Permut,j]=find(P');
end