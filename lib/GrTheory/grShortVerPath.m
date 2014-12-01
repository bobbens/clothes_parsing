function [dMWP,ssp]=grShortVerPath(E,Wv)
% Function dMVP=grShortVerPath(E,Wv) for digraph with weighted vertexes 
% solve the problem about the path with minimal weight of verticies.
% Input parameters: 
%   E(m,2) - the arrows of digraph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of arrows.
%   Wv(n,1) - the weights of verticies; if this parameter omitted, 
%             then all vertices have weight 1.
%     n - number of verticies.
% Output parameter:
%   dMWP(k) - the vector with number of verticies included
%     to path with minimal weight.
% [dMWP,ssp]=grMaxVerPath(E,Wv) return also
%   ssp - sum of vertexes weights on this path.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru
% Acknowledgements to Dr. Albert Niel (Austria)
% for testing of this algorithm.

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation
if nargin<2,
  Wv = ones(n,1);
end
if n~=length(Wv),
  error('Length of Wv not equal to n!')
end
%============ Duplication of vertices ==================
E1 = [[E(:,1)+n, E(:,2), zeros(m,1)]; [(1:n)', (n+1:2*n)', Wv(:)]];
%============ Finding of the base and contrabase =======
BG = grBase(E1); % base of digraph
CBG = grCoBase(E1); % contrabase of digraph
S1 = unique(BG(:))'; % all vertices-sources
T1 = unique(CBG(:))'; % all vertices-tails
ns = length(S1); % number of vertices in base
nt = length(T1); % number of vertices in contrabase
%============ Add source and tail vertexes =============
s = 2*n+1; % source vertex
t = 2*n+2; % tail vertex
E2=[[ones(1,ns)*s;S1;zeros(1,ns)]';E1;[T1;ones(1,nt)*t;zeros(1,nt)]'];
%============ Shortest Path Problem ====================
[dSP,sp] = grShortPath(E2,s,t);
dMWP = sp(2:2:end-1);
ssp = sum(Wv(dMWP));
return