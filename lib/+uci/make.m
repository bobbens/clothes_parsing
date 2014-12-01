function make(varargin)
%MAKE build uci package

target = 'all';
if nargin > 0, target = varargin{1}; end

cwd = pwd;
ROOT = fullfile(fileparts(mfilename('fullpath')),'private');
cd(ROOT);
try
switch target
    case 'all'
        mex -O all_pairs_bb_iou.cc
        mex -O vgg_nearest_neighbour_dist.cc
        mex -O qp_one_sparse.cc
        mex -O score.cc
        mex -O resize.cc
        mex -O reduce.cc
        mex -O dt.cc
        mex -O features.cc
        mex -O fconv.cc
    case 'clean'
        delete *.mex*
end
catch e
    disp(e.getReport);
end
cd(cwd);

end

