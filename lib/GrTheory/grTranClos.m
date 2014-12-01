function Etc=grTranClos(E)
% Function Etc=grTranClos(E) built 
% the transitive closure for the digraph E.
% Input parameter: 
%   E(m,2) - the arrows of digraph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of arrows.
% Output parameter:
%   Etc(mtc,2) - the arrows of digraph with trasitive closure for E;
%     mtc - number of arrows.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % data validation
% ============ Atc=A+A^2+A^3+A^4+... ==========================
A1=zeros(n);
A1((E(:,2)-1)*n+E(:,1))=1; % adjacency matrix A
A0=A1;
An=(A1*A1)>0; % A^2
A=(A1+An)>0; % A+A^2
while ~all(all(A==A1)),
  A1=A;
  An=((double(An))*A0)>0;
  A=(A+An)>0;
end
[r,c]=find(A);
Etc=[r c]; % all arrows of the transitive closure
return