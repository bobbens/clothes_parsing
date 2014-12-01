%MATLAB_INIT    Initialize Matlab environment.

% add paths
p = fileparts(mfilename('fullpath'));
path(path,fullfile(p,''));
path(path,fullfile(p,'util',''));
path(path,fullfile(p,'util','math',''));
path(path,fullfile(p,'util','ui',''));
