function [Ec,Rad,Diam,Cv,Pv]=grEccentricity(E)
% Function Ec=grEccentricity(E) find the (weighted) 
% eccentricity of all vertexes of graph.
% Input parameter: 
%   E(m,2) or (m,3) - the edges of graph and their weight;
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of arrow;
%     m - number of arrows.
%     If we set the array E(m,2), then all weights is 1.
% Output parameter:
%   Ec(1,n) - the (weighted) eccentricity of all vertexes.
% [Ec,Rad,Diam]=grEccentricity(E) find also the radius Rad and 
%   diameter Diam of the graph (the scalars).
% [Ec,Rad,Diam,Cv,Pv]=grEccentricity(E) find also 
%   the center vertexes Cv and the periphery vertexes Pv 
%   of the graph (the vector-rows with numbers of vertexes).
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

if nargin<1,
  error('There are no input data!')
end
dSP=grDistances(E); % the matrix of distances
Ec=max(dSP); % the eccentricity of all vertexes
if nargout>=2,
  Rad=min(Ec); % the radius
  if nargout>=3,
    Diam=max(Ec); % the diameter
    if nargout>=4,
      Cv=find(Ec==Rad); % the center vertexes of the graph
      if nargout>=5,
        Pv=find(Ec==Diam); % the periphery vertexes of the graph
      end
    end
  end
end
return