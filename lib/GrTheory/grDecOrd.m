function [Dec,Ord]=grDecOrd(E)
% Function Dec=grDecOrd(E) solve 
% the problem about decomposition of the digraph
% to the sections with mutually accessed vertexes
% (strongly connected components).
% Input parameter: 
%   E(m,2) - the arrows of digraph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of arrows.
% Output parameter:
%   Dec(n,ns) - the Boolean array with numbers of vertexes.
%     n - number of vertexes;
%     ns - number of sections with mutually accessed vertexes.
%     In each column of the array Dec True value have
%     numbers of vertexes of this section.
% Other syntax: [Dec,Ord]=grDecOrd(E) also ordered all sections
%   in partial ordering. Second output parameter:
%   Ord(ns,ns) - the Boolean matrix of partial ordering
%     of sections. This matrix have right-up triangle structure.
%     If Ord(i,j)=True, then we have access 
%     from the section i (vertexes of i-st column of Dec)
%     to the section j (vertexes of j-st column of Dec).
%     In this syntax all columns of Dec(n,ns) 
%     is in partial ordering.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru
% Acknowledgements to Dr.Eduardo D. Sontag (sontag@math.rutgers.edu)
% for testing of this algorithm.

if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % data validation
% ============ Decomposition ==========================
A=eye(n);
A((E(:,2)-1)*n+E(:,1))=1; % adjacency matrix with main diagonal
A2=(A*A)>0;
while ~all(all(A2==A)),
  A=double(A2);
  A2=(A*A)>0;
end
[T,ir,jc]=unique(A.*A','rows');
Dec=T';

% ============== Partial ordering =========================
if nargout>1, 
  Ord=A(ir,ir);
  ns=size(Ord,1); % number of sections
  for it=1:ns*(ns-1)/2, % the iterations for partial ordering
    Mlow=tril(Ord,-1); % the left down trialgle
    [is,js,Mw]=find(Mlow); % we find not ordering elements
    if isempty(is), % all ordered
      break; % exit from loop for
    end
    num=[1:ns];
    num(is(1))=js(1); % we change two numbers
    num(js(1))=is(1);
    Ord=Ord(num,num);
    Dec=Dec(:,num);
  end
end  
return