function model = trainmodel(name,pos,neg,K,pa,sbin)

%globals;

cls = [name '_cluster_' num2str(K')' '.mat'];
if exist(fullfile(cachedir,cls),'file')
  load(fullfile(cachedir,cls));
else
  model = initmodel(pos,sbin);
  def = data_def(pos,model);
  idx = clusterparts(def,K,pa);
  save(fullfile(cachedir,cls),'def','idx');
end

numpart = length(pa);
for p = 1:numpart
  cls = [name '_part_' num2str(p) '_mix_' num2str(K(p)) '.mat'];
  if exist(fullfile(cachedir,cls),'file')
    load(fullfile(cachedir,cls));
  else
    sneg = neg(1:min(length(neg),100));
    model = initmodel(pos,sbin);
    models = cell(1,K(p));
    for k = 1:K(p)
      spos = pos(idx{p} == k);
      for n = 1:length(spos)
        spos(n).x1 = spos(n).x1(p);
        spos(n).y1 = spos(n).y1(p);
        spos(n).x2 = spos(n).x2(p);
        spos(n).y2 = spos(n).y2(p);
      end
      models{k} = train(cls,model,spos,spos,1,1);
      % models{k} = train(cls,model,spos,sneg,1,1);
    end
    model = mergemodels(models);
    save(fullfile(cachedir,cls),'model');
  end
end


cls = [name '_final_' num2str(K')' '.mat'];
if exist(fullfile(cachedir,cls),'file')
  load(fullfile(cachedir,cls));
else
  model = buildmodel(name,model,pa,def,idx,K);
  %model = train(cls,model,pos,neg,0,2);
  model = train(cls,model,pos,pos,0,2);
  save(fullfile(cachedir,cls),'model');
end
end

%%
function deffeat = data_def(pos,model)
width  = zeros(1,length(pos));
height = zeros(1,length(pos));
labels = zeros(size(pos(1).point,1),size(pos(1).point,2),length(pos));
for n = 1:length(pos)
  width(n)  = pos(n).x2(1) - pos(n).x1(1);
  height(n) = pos(n).y2(1) - pos(n).y1(1);
  labels(:,:,n) = pos(n).point;
end
scale = sqrt(width.*height)/sqrt(model.maxsize(1)*model.maxsize(2));
scale = [scale; scale];

deffeat = cell(1,size(labels,1));
def0 = squeeze(labels(1,1:2,:));
for p = 1:size(labels,1)
  def = squeeze(labels(p,1:2,:));
  def = (def - def0) ./ scale;
  deffeat{p} = def';
end
end

%%
function jointmodel = buildmodel(name,model,pa,def,idx,K)
% jointmodel = buildmodel(name,model,pa,def,idx,K)
% This function merges together separate part models into a tree structure

%globals;

jointmodel.bias    = struct('w',{},'i',{});
jointmodel.defs    = struct('w',{},'i',{},'anchor',{});
jointmodel.filters = struct('w',{},'i',{});
jointmodel.components{1} = struct('biasid',{},'defid',{},'filterid',{},'parent',{});

jointmodel.maxsize  = model.maxsize;
jointmodel.interval = model.interval;
jointmodel.sbin = model.sbin;
jointmodel.len = 0;

% add children
for i = 1:length(pa)
  child = i;
  parent = pa(child);
  assert(parent < child);
  
  cls = [name '_part_' num2str(child) '_mix_' num2str(K(i)) '.mat'];
  load(fullfile(cachedir,cls));

  % add bias
  p.biasid = [];
  if parent == 0
    nb  = length(jointmodel.bias);
    b.w = 0;
    b.i = jointmodel.len + 1;
    jointmodel.bias(nb+1) = b;
    jointmodel.len = jointmodel.len + numel(b.w);
    p.biasid = nb+1;
  else
    for k = 1:max(idx{child})
      for l = 1:max(idx{parent})
        % if any(idx{child} == k & idx{parent} == l),
        nb = length(jointmodel.bias);
        b.w = 0;
        b.i = jointmodel.len + 1;
        jointmodel.bias(nb+1) = b;
        jointmodel.len = jointmodel.len + numel(b.w);
        p.biasid(l,k) = nb+1;
      end
    end
  end

  % add deformation parameter
  p.defid = [];
  if parent > 0
    for k = 1:max(idx{child})
      nd  = length(jointmodel.defs);
      d.w = [0.01 0 0.01 0];
      d.i = jointmodel.len + 1;
      x = mean(def{child}(idx{child}==k,1) - def{parent}(idx{child}==k,1)); 
      y = mean(def{child}(idx{child}==k,2) - def{parent}(idx{child}==k,2));
      d.anchor = round([x+1 y+1 0]);
      jointmodel.defs(nd+1) = d;
      jointmodel.len = jointmodel.len + numel(d.w);	
      p.defid = [p.defid nd+1];
    end
  end

  % add filter
  p.filterid = [];
  for k = 1:max(idx{child})    
    nf  = length(jointmodel.filters);
    f.w = model.filters(k).w;
    f.i = jointmodel.len + 1;
    jointmodel.filters(nf+1) = f;
    jointmodel.len = jointmodel.len + numel(f.w);
    p.filterid = [p.filterid nf+1];
  end

  p.parent = parent;
  np = length(jointmodel.components{1});
  jointmodel.components{1}(np+1) = p;
end
end