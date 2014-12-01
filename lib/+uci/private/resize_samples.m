function [ X, s ] = resize_samples( X, siz )
%RESIZE_SAMPLES resize samples to a specified dimension

s = ones(size(X));
for i = 1:numel(X)
    if ischar(X(i).im), X(i).im = imread(X(i).im); end
    d = size(X(i).im);
    s(i) = min(siz./d(1:2)); % scaling factor
    if s(i) < 1
        X(i).im = imresize(X(i).im,s(i));
        if isfield(X(i),'point')
            X(i).point = (X(i).point-1).*s(i) + 1;
        end
        if ~isempty(X(i).labels)
            X(i).labels = imresize(X(i).labels,s(i),'nearest');
        end
    end
end
assert(all(s>0));

end

