function pyra = featpyramid(S, model)
% Compute feature pyramid.
%
% pyra.feat{i} is the i-th level of the feature pyramid.
% pyra.scales{i} is the scaling factor used for the i-th level.
% pyra.feat{i+interval} is computed at exactly half the resolution of feat{i}.
% first octave halucinates higher resolution data.

% BEGIN HACK %%%%%%%%%%%%%%%
im = S.im;
map = S.labels;
if ischar(im), im = imread(im); end
% END HACK %%%%%%%%%%%%%%%%

sbin      = model.sbin;
interval  = model.interval;
padx      = max(model.maxsize(2)-1-1,0);
pady      = max(model.maxsize(1)-1-1,0);
sc = 2 ^(1/interval);
imsize = [size(im, 1) size(im, 2)];
max_scale = 1 + floor(log(min(imsize)/(5*sbin))/log(sc));
pyra.feat = cell(max_scale,1);
pyra.scale = zeros(max_scale,1);

if size(im, 3) == 1
  im = repmat(im,[1 1 3]);
end
im = double(im); % our resize function wants floating point values

for i = 1:interval
  scaled = resize(im, 1/sc^(i-1));
  pyra.feat{i} = features(scaled,sbin);
  % BEGIN HACK %%%%%%%%%%%%%%%%%%
  nsift = size(pyra.feat{i},3);
  if ~isempty(map)
      map = imresize(map,[size(scaled,1) size(scaled,2)],'nearest');
      pyra.feat{i} = cat(3,pyra.feat{i},build_cloth_features(map,sbin));
  end
  % END HACK %%%%%%%%%%%%%%%%%%%%
  pyra.scale(i) = 1/sc^(i-1);
  % remaining interals
  for j = i+interval:interval:max_scale
    scaled = reduce(scaled);
    pyra.feat{j} = features(scaled,sbin);
    % BEGIN HACK %%%%%%%%%%%%%%%%%%
    if ~isempty(map)
      map = imresize(map,[size(scaled,1) size(scaled,2)],'nearest');
      pyra.feat{j} = cat(3,pyra.feat{j}, build_cloth_features(map,sbin));
    end
    % END HACK %%%%%%%%%%%%%%%%%%%%
    pyra.scale(j) = 0.5 * pyra.scale(j-interval);
  end
end

for i = 1:length(pyra.feat)
  % add 1 to padding because feature generation deletes a 1-cell
  % wide border around the feature map
  pyra.feat{i} = padarray(pyra.feat{i}, [pady+1 padx+1 0], 0);
  % write boundary occlusion feature
  pyra.feat{i}(1:pady+1, :, nsift) = 1;
  pyra.feat{i}(end-pady:end, :, nsift) = 1;
  pyra.feat{i}(:, 1:padx+1, nsift) = 1;
  pyra.feat{i}(:, end-padx:end, nsift) = 1;
end

pyra.scale    = model.sbin./pyra.scale;
pyra.interval = interval;
pyra.imy = imsize(1);
pyra.imx = imsize(2);
pyra.pady = pady;
pyra.padx = padx;

end
