function [warped lwarped] = warppos(name, model, pos, verbose)
% warped = warppos(name, model, pos)
% Warp positive examples to fit model dimensions.
% Used for training root filters from positive bounding boxes.

if nargin < 4, verbose = true; end

f   = model.components{1}(1).filterid;
siz = size(model.filters(f).w);
siz = siz(1:2);
pixels = siz * model.sbin; 
heights = [pos(:).y2]' - [pos(:).y1]' + 1;
widths = [pos(:).x2]' - [pos(:).x1]' + 1;
numpos = length(pos);
warped = cell(numpos,1);
% BEGIN HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lwarped = cell(numpos,1);
% END HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cropsize = (siz+2) * model.sbin;
for i = 1:numpos
%   fprintf('%s: warp: %d/%d\n', name, i, numpos);
%   im = imread(pos(i).im);
%   if size(im, 3) == 1
%     im = repmat(im,[1 1 3]);
%   end
  % BEGIN HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if verbose, fprintf('%s: warp: %d/%d\n', name, i, numpos); end
  im = pos(i).im;
  if ischar(im), im = imread(im); end
  % END HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  padx = model.sbin * widths(i) / pixels(2);
  pady = model.sbin * heights(i) / pixels(1);
  x1 = round(pos(i).x1-padx);
  x2 = round(pos(i).x2+padx);
  y1 = round(pos(i).y1-pady);
  y2 = round(pos(i).y2+pady);
  window = subarray(im, y1, y2, x1, x2, 1);
  warped{i} = imresize(window, cropsize, 'bilinear');
  % BEGIN HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~isempty(pos(i).labels)
    window = subarray(pos(i).labels,y1,y2,x1,x2,0);
    lwarped{i} = imresize(window,cropsize,'nearest');
  end
  % END HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function B = subarray(A, i1, i2, j1, j2, pad)

% B = subarray(A, i1, i2, j1, j2, pad)
% Extract subarray from array
% pad with boundary values if pad = 1
% pad with zeros if pad = 0

dim = size(A);
if numel(dim)==2, dim = [dim,1]; end
B = zeros(i2-i1+1, j2-j1+1, dim(3));
if pad
  for i = i1:i2
    for j = j1:j2
      ii = min(max(i, 1), dim(1));
      jj = min(max(j, 1), dim(2));
      B(i-i1+1, j-j1+1, :) = A(ii, jj, :);
    end
  end
else
  for i = max(i1,1):min(i2,dim(1))
    for j = max(j1,1):min(j2,dim(2))
      B(i-i1+1, j-j1+1, :) = A(i, j, :);
    end
  end
end

