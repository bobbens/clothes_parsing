function model = train(name, model, pos, neg, warp, iter, C, wpos, maxsize, overlap, verbose) 
% model = train(name, model, pos, neg, warp, iter, C, Jpos, maxsize, overlap)
%               1,    2,     3,   4,   5,    6,    7, 8,    9,       10
% Train a structured SVM with latent assignement of positive variables
% pos  = list of positive images with part annotations
% neg  = list of negative images
% warp =1 uses warped positives
% warp =0 uses latent positives
% iter is the number of training iterations
%   C  = scale factor for slack loss
% wpos =  amount to weight errors on positives
% maxsize = maximum size of the training data cache (in GB)
% overlap =  minimum overlap in latent positive search
if nargin < 11, verbose = true; end

if nargin < 6
  iter = 1;
end

if nargin < 7
  C = 0.002;
end

if nargin < 8
  wpos = 2;
end

if nargin < 9
  % Estimated #sv = (wpos + 1) * # of positive examples
  % maxsize*1e9/(4*model.len)  = # of examples we can store, encoded as 4-byte floats
  no_sv = (wpos+1) * length(pos);
  maxsize = 10 * no_sv * 4 * sparselen(model) / 1e9;
  maxsize = min(max(maxsize,3),6);
  %maxsize = 3; %Uncomment this line to run comfortably on machines with 4GB of memory
  %maxsize = 5; %Uncomment this line to run comfortable on machines with 8GB memory
end

if verbose, fprintf('Using %.1f GB\n',maxsize); end

if nargin < 10
  overlap = 0.6;
end


% Vectorize the model
len  = sparselen(model);
nmax = round(maxsize*.25e9/len);

%rand('state',0);
%globals;

% Define global QP problem
clear global qp;
global qp;
% qp.x(:,i) = examples
% qp.i(:,i) = id
% qp.b(:,i) = bias of linear constraint
% qp.d(i)   = ||qp.x(:,i)||^2
% qp.a(i)   = ith dual variable
qp.x   = zeros(len,nmax,'single');
qp.i   = zeros(5,nmax,'int32');
qp.b   = ones(nmax,1,'single');
qp.d   = zeros(nmax,1,'double');
qp.a   = zeros(nmax,1,'double');
qp.sv  = false(1,nmax);  
qp.n   = 0;
qp.obj = [];

[qp.w,qp.wreg,qp.w0,qp.noneg] = model2vec(model);
qp.Cpos = C*wpos;
qp.Cneg = C;
qp.w    = (qp.w - qp.w0).*qp.wreg;

for t = 1:iter,
  if verbose, fprintf('\niter: %d/%d\n', t, iter); end
  qp.n = 0;
  if warp > 0
    numpositives = poswarp(name, t, model, pos);
  else
    numpositives = poslatent(name, t, model, pos, overlap);
  end
  
  if verbose
      for i = 1:length(numpositives)
        fprintf('component %d got %d positives\n', i, numpositives(i));
      end
  end
  assert(qp.n <= nmax);
  
  % Fix positive examples as permenant support vectors
  % Initialize QP soln to a valid weight vector
  % Update QP with coordinate descent
  qp.svfix = 1:qp.n;
  qp.sv(qp.svfix) = 1;
  qp_prune();
  qp_opt(1);
  model = vec2model(qp_w,model);
  model.interval = 2;
  loss = 0;
  
  % grab negative examples from negative images
  for i = 1:length(neg),
    if verbose, fprintf('\n Image(%d/%d)',i,length(neg)); end
    % BEGIN HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bbox = [];
    if isfield(neg(i),'point')
        bbox = [neg(i).x1;neg(i).y1;neg(i).x2;neg(i).y2];
    end
    [box,model,loss] = detect(neg(i), model, -1, bbox(:)', overlap, i, -1,loss);
    % END HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if verbose, fprintf(' #cache+%d=%d/%d, #sv=%d, LB=%.4f, loss=%.4f',...
            size(box,1),qp.n,nmax,sum(qp.sv),qp.obj,loss); end
    % Stop if cache is full
    if sum(qp.sv) == nmax, break; end
  end

  % One final pass of optimization
  qp_opt();
  model = vec2model(qp_w,model);

  if verbose, fprintf('\nDONE iter: %d/%d #sv=%d/%d, LB=%.4f\n',...
          t,iter,sum(qp.sv),nmax,qp.obj); end

  % Compute minimum score on positive example (with raw, unscaled features)
  r = sort(qp_scorepos);
  model.thresh   = r(ceil(length(r)*.05));
  model.interval = 10;
  model.obj = qp.obj;
  % visualizemodel(model);
  % cache model
  % save([cachedir name '_model_' num2str(t)], 'model');
end
if verbose, fprintf('qp.x size = [%d %d]\n',size(qp.x)); end
clear global qp;
end

% get positive examples by warping positive bounding boxes
% we create virtual examples by flipping each image left to right
function numpositives = poswarp(name, t, model, pos, verbose)
  if nargin < 5, verbose = true; end
  numpos = length(pos);
  %BEGIN HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [warped, lwarped] = warppos(name, model, pos);
  % END HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  minsize = prod(model.maxsize*model.sbin);

  for i = 1:numpos
    if verbose, fprintf('%s: iter %d: warped positive: %d/%d\n',...
            name, t, i, numpos); end
    bbox = [pos(i).x1 pos(i).y1 pos(i).x2 pos(i).y2];
    % skip small examples
    if (bbox(3)-bbox(1)+1)*(bbox(4)-bbox(2)+1) < minsize
      continue
    end    
    % get example
    im = warped{i};
    feat = features(im, model.sbin);
    % BEGIN HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    feat = cat(3,feat,build_cloth_features(lwarped{i},model.sbin));
    % END HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qp_poswrite(feat,i,model);
  end
  global qp;
  numpositives = qp.n;
end
  
function qp_poswrite(feat,id,model)

  len = numel(feat);
  ex.id     = [1 id 0 0 0]';
  ex.blocks = [];
  ex.blocks(end+1).i = model.bias.i;
  ex.blocks(end).x   = 1;
  ex.blocks(end+1).i = model.filters.i;
  ex.blocks(end).x   = feat;
  qp_write(ex);
end

% get positive examples using latent detections
% we create virtual examples by flipping each image left to right
function numpositives = poslatent(name, t, model, pos, overlap, verbose)
  if nargin < 6, verbose = true; end
  numpos = length(pos);
  model.interval = 5;
  numpositives = zeros(length(model.components), 1);
  minsize = prod(model.maxsize*model.sbin);
  
  % we add the below lines 
  numparts = length(model.components{1});
  for i = 1:numpos
    if verbose, fprintf('%s: iter %d: latent positive: %d/%d\n',...
            name, t, i, numpos); end
    % skip small examples
    skipflag = 0;
    bbox = cell(1,numparts);
    for p = 1:numparts
      bbox{p} = [pos(i).x1(p) pos(i).y1(p) pos(i).x2(p) pos(i).y2(p)];
      if (bbox{p}(3)-bbox{p}(1)+1)*(bbox{p}(4)-bbox{p}(2)+1) < minsize
        skipflag = 1;
        break;
      end
    end
    if skipflag
      continue;
    end
    
    % get example
%     im = imread(pos(i).im);
%     [im, bbox] = croppos(im, bbox);
    % BEGIN HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    im = pos(i).im;
    if ischar(im), im = imread(pos(i).im); end
    [im, bbox] = croppos(im, bbox);
    cropped_labels = croppos(pos(i).labels, bbox);
    X = struct('im',im,'labels',cropped_labels);
    % END HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    box = detect(X, model, 0, bbox, overlap, i, 1);
    if ~isempty(box)
      if verbose, fprintf(' (comp=%d,sc=%.3f)\n',...
              box(1,end-1),box(1,end)); end
      c = box(1,end-1);
      numpositives(c) = numpositives(c)+1;
      %      showboxes(im, box);
    end
  end
end

% Compute score (weights*x) on positives examples (see qp_write.m)
% Standardized QP stores w*x' where w = (weights-w0)*r, x' = c_i*(x/r)
% (w/r + w0)*(x'*r/c_i) = (v + w0*r)*x'/ C
function scores = qp_scorepos
  global qp;
  y = qp.i(1,1:qp.n);
  I = find(y == 1);
  w = qp.w + qp.w0.*qp.wreg;
  scores = score(w,qp.x,I) / qp.Cpos;
end

% Computes expected number of nonzeros in sparse feature vector 
function len = sparselen(model)

  len       = 0;
  numblocks = 0;
  for c = 1:length(model.components),
    feat = zeros(model.len,1);
    for p = model.components{c},
      if ~isempty(p.biasid)
        x = model.bias(p.biasid(1));
        i1 = x.i;
        i2 = i1 + numel(x.w) - 1;
        feat(i1:i2) = 1;
        numblocks = numblocks + 1;
      end
      if ~isempty(p.defid)
        x  = model.defs(p.defid(1));
        i1 = x.i;
        i2 = i1 + numel(x.w) - 1;
        feat(i1:i2) = 1;
        numblocks = numblocks + 1;
      end
      if ~isempty(p.filterid)
        x  = model.filters(p.filterid(1));
        i1 = x.i;
        i2 = i1 + numel(x.w) - 1;
        feat(i1:i2) = 1;
        numblocks = numblocks + 1;
      end
    end
    
    % Number of entries needed to encode a block-sparse representation
    %   1 + numberofblocks*2 + #nonzeronumbers
    n = 1 + numblocks*2 + sum(feat);
    len = max(len,n);
  end
end

% [newim, newbox] = croppos(im, box)
% Crop positive example to speed up latent search.
function [im, box] = croppos(im, box)

  if iscell(box) % many parts are stored in box
    P = length(box);
    x1 = zeros(1,P);
    y1 = zeros(1,P);
    x2 = zeros(1,P);
    y2 = zeros(1,P);
    for p = 1:P
      x1(p) = box{p}(1);
      y1(p) = box{p}(2);
      x2(p) = box{p}(3);
      y2(p) = box{p}(4);
    end
    x1 = min(x1); y1 = min(y1); x2 = max(x2); y2 = max(y2);
    pad = 0.5*((x2-x1+1)+(y2-y1+1));
    x1 = max(1, round(x1-pad));
    y1 = max(1, round(y1-pad));
    x2 = min(size(im,2), round(x2+pad));
    y2 = min(size(im,1), round(y2+pad));
    
    im = im(y1:y2, x1:x2, :);
    for p = 1:P
      box{p}([1 3]) = box{p}([1 3]) - x1 + 1;
      box{p}([2 4]) = box{p}([2 4]) - y1 + 1;
    end
    
  else % only one template is stored in box
       % crop image around bounding box
    pad = 0.5*((box(3)-box(1)+1)+(box(4)-box(2)+1));
    x1 = max(1, round(box(1) - pad));
    y1 = max(1, round(box(2) - pad));
    x2 = min(size(im, 2), round(box(3) + pad));
    y2 = min(size(im, 1), round(box(4) + pad));
    
    im = im(y1:y2, x1:x2, :);
    box([1 3]) = box([1 3]) - x1 + 1;
    box([2 4]) = box([2 4]) - y1 + 1;
  end
end
