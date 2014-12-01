function [detRate PCP R] = eval_pcp(name,predstickmen,gtstickmen,thresh)
% Please refer to BatchEvalBuffy.m in buffy_stickmen_v2.1/code/ for
% understanding

if nargin < 4
  thresh  = 0.5;
end

nLimbs = size(gtstickmen(1).stickmen.coor,2);

nTotal = length(predstickmen);

assert(nTotal == length(gtstickmen));

scores = zeros(1,nTotal);
Rs = zeros(nTotal,nLimbs);

for i = 1:nTotal
  [trash scores(i) rs] = EvalStickmen(name,predstickmen(i).stickmen,gtstickmen(i).stickmen,thresh);
  if ~isempty(rs)
    Rs(i,:) = rs;
  end
end

% compute detRate and PCP
matched = ~isnan(scores);
nMatched = sum(matched);
detRate = nMatched/nTotal;
PCP = sum(scores(matched))/(nMatched*nLimbs);
R = sum(Rs(matched,:))/(nMatched);
end