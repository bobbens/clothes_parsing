function nCol=grColVer(E)
% function nCol=grColVer(E) solve the color graph problem
% for vertexes of the graph.
% Input parameter: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
% Output parameter:
%   nCol(n,1) - the list of the colors of vertexes.
% Uses the reduction to integer LP-problem.
% Required the Optimization Toolbox v.3.0.1 or over.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation

E=sort(E(:,1:2)')'; % each row in ascending order
E=unique(E,'rows'); % we delete multiple edges
E=E(setdiff([1:size(E,1)]',find((E(:,1)==E(:,2)))),:); % we delete loops

% ============= Parameters of integer LP problem ==========
A1=zeros(2*m,n*n);
for kk=1:m,
  A1(2*kk-1,(E(kk,1)-1)*n+1:E(kk,1)*n)=1;
  A1(2*kk-1,(E(kk,2)-1)*n+1:E(kk,2)*n)=-1;
  A1(2*kk,(E(kk,1)-1)*n+1:E(kk,1)*n)=-1;
  A1(2*kk,(E(kk,2)-1)*n+1:E(kk,2)*n)=1;
end
A=[zeros(2*m,n),A1,...
  reshape([-n*reshape(eye(m),1,m*m);n*reshape(eye(m),1,m*m)],2*m,m);...
  -ones(n),reshape(repmat(reshape(eye(n),1,n*n),n,1),n*n,n)',zeros(n,m)];
b=[reshape([-ones(1,m);(n-1)*ones(1,m)],2*m,1);zeros(n,1)];
c=[ones(n,1);zeros(n*n+m,1)];
options=optimset('bintprog'); % the default options
options.Display='off'; % we change the output

% ============= We solve the integer LP problem ==========
xmin=round(bintprog(c,A,b,[],[],[],options));
nCol=(sum(reshape(xmin(n+1:n*(n+1)),n,n))+1)'
return