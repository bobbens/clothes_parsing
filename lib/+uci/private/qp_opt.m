function qp_opt(tol)
% qp_opt(tol)
% Optimize QP until relative difference between lower and upper bound is below 'tol'
% If tol = 1, then optimization will stop when lower bound stops increasing

global qp;

if nargin < 1,
  tol = .05;
end

C = 1;
I = 1:qp.n;
[id,J] = sortrows(qp.i(:,I)');
id     = id';
eqid   = [0 all(id(:,2:end) == id(:,1:end-1),1)];

slack = qp.b(I) - score(qp.w,qp.x,I);
loss  = computeloss(slack(J),eqid);
ub    = qp.w'*qp.w*.5 + C*loss; 
lb    = qp.obj;
fprintf('\n LB=%.4f,UB=%.4f [',lb,ub);
for t = 1:100,
  % To avoid computing upper bounds,
  % iterate until there is no improvement on the dual 
  % when starting with a full active set
  for itr = 1:10 % limit iteration to 10
    i   = 0;
    obj = -inf;
    qp.sv(1:qp.n) = 1;
    while qp.obj < 0 || ((qp.obj - obj)/qp.obj > .001)
      obj = qp.obj;
      qp_one;
      fprintf('.');
      i = i + 1;
      if i >= 10, break; end
    end
    if i <= 1, break; end
  end
  if itr >= 100, warning('qp_opt:maxitr',...
      'Number of qp_one iteration reached 10'); end
  
  % Compute upper and lower bounds
  slack = qp.b(I) - score(qp.w,qp.x,I);
  loss  = computeloss(slack(J),eqid);
  ub    = min(ub,qp.w'*qp.w*.5 + C*loss);  
  lb    = qp.obj;  
  if 1 - lb/ub < tol,
    break;
  end
  %fprintf('t=%d: UB=%.5f,LB=%.5f\n',t,ub,lb);
end

fprintf('] LB=%.4f,UB=%.4f\n',lb,ub);
end


function loss = computeloss(slack,eqid)
% Zero-out scores that aren't the greatest violated constraint for an id
% eqid(i) = 1 if x(i) and x(i-1) are from the same id
% eqid(1) = 0
% v is the best value in the current block
% i is a ptr to v
% j is a ptr to the example we are considering

err = false(size(eqid));
for j = 1:length(err),
  % Are we at a new id?
  % If so, update i,v
  if eqid(j) == 0,
    i = j;
    v = slack(i);
    if v > 0,
      err(i) = 1;
    else
      v = 0;
    end
    % Are we at a new best in this set of ids?
    % If so, update i,v and zero out previous best
  elseif slack(j) > v 
    err(i) = 0;
    i = j;
    v = slack(i);
    err(i) = 1;
  end
end

loss = sum(slack(err));
end
