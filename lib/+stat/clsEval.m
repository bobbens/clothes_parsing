function [S, order] = clsEval( X, Xhat, varargin )
%CLSEVAL calculate evaluation metrics for a multi-way classification
%
%    [S, order] = stat.clsEval(X, Xhat)
%    [...] = stat.clsEval(..., 'param1', value1, ...)
%
% # Input
% * __X__ true label vector
% * __Xhat__ predicted label vector
% 
% # Params
% * __Order__ order of unique labels in the confusion matrix. By default,
%   it is set to union(X,Xhat).
% * __SampleWeights__ weights applied to samples
% * __ClassWeights__ weights applied to classes
%

% Inspect input
X = X(:);
Xhat = Xhat(:);
assert(numel(X)==numel(Xhat));

% Set options
W = ones(size(X));
order = [];
CW = [];
ignore_nan = true;
for i = 1:2:numel(varargin)
    switch varargin{i}
        case 'Order', order = varargin{i+1};
        case 'SampleWeights', W = varargin{i+1};
        case 'ClassWeights', CW = varargin{i+1};
        case 'IgnoreNaN', ignore_nan = varargin{i+1};
    end
end
if isempty(order), order = union(X,Xhat); end
if isempty(CW), CW = ones(size(order)); end
W = W(:);
CW = CW(:);
sumW = sum(W);
sumCW = sum(CW);

assert(numel(CW)==numel(order));
assert(numel(W)==numel(X));
assert(sumW~=0);
assert(sumCW~=0);

% Index labels
if iscellstr(X) && iscellstr(Xhat)
    Y = arrayfun(@(x) find(strcmp(order,x),1), X);
    Yhat = arrayfun(@(x) find(strcmp(order,x),1), Xhat);
else
    Y = arrayfun(@(x) find(order==x,1), X);
    Yhat = arrayfun(@(x) find(order==x,1), Xhat);
end

% Calculate evaluations
S.acc = sum(W.*(Y(:)==Yhat(:))) / sumW;
S.confusion = full(sparse(Y,Yhat,W,numel(order),numel(order)));

% Additional measures
S.pre = diag(S.confusion) ./ sum(S.confusion,1)';
S.rec = diag(S.confusion) ./ sum(S.confusion,2);
S.f1  = 2*(S.pre.*S.rec)   ./ (S.pre+S.rec);
if ignore_nan
    ind = ~isnan(S.pre);
    S.avr_pre = sum(CW(ind).*S.pre(ind)) / sum(CW(ind));
    ind = ~isnan(S.rec);
    S.avr_rec = sum(CW(ind).*S.rec(ind)) / sum(CW(ind));
    ind = ~isnan(S.f1);
    S.avr_f1  = sum(CW(ind).*S.f1(ind))  / sum(CW(ind));
else
    S.avr_pre = sum(CW.*S.pre) / sumCW;
    S.avr_rec = sum(CW.*S.rec) / sumCW;
    S.avr_f1  = sum(CW.*S.f1)  / sumCW;
end
S.order = order;

end

