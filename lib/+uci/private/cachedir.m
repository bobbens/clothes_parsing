function r = cachedir(q)
%CACHEDIR path to cachedir
persistent p;
if isempty(p), p = fullfile('cache'); end
if nargin > 0, p = q; end

if ~exist(p,'dir'), mkdir(p); end

r = p;
end