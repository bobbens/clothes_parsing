% load exemplar data from the given file
function [features, x_pos, y_pos] = load_exemplar(filename, ftype)

% check whether to use pca features
if (strcmp(ftype,'no_pca'))
   use_pca = 0;
elseif (strcmp(ftype,'pca'))
   use_pca = 1;
else
   error('ftype must be no_pca or pca');
end

% load features
data = load(filename);
if (use_pca)
   features = data.pv./repmat(sqrt(sum(data.pv .* data.pv,2)),[1 size(data.pv,2)]);
else
   n_desc = size(data.descriptors,1);
   n_dim = size(data.descriptors,2)*size(data.descriptors,3);
   features = reshape(data.descriptors,[n_desc n_dim]);
end
x_pos = data.pos(:,1);
y_pos = data.pos(:,2);
