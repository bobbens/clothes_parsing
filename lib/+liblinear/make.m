function make(target)
%MAKE build liblinear matlab package
%
%  liblinear.make
%  liblinear.make('all')
%  liblinear.make('clean')
%
    if nargin < 1, target = 'all'; end
    cwd = pwd;
    p = fileparts(mfilename('fullpath'));
    cd(p);
    switch computer('arch')
        case {'glnx86','glnxa64','maci','maci64'}
            system(['make MATLABDIR=',matlabroot,' ',target]);
        case {'win32','win64'}
            % add -largeArrayDims on 64-bit machines
            mex -O -largeArrayDims -c blas\*.c -outdir blas
            mex -O -largeArrayDims -c linear.cpp
            mex -O -largeArrayDims -c tron.cpp
            mex -O -largeArrayDims -c linear_model_matlab.c
            mex -O -largeArrayDims train.c tron.obj linear.obj linear_model_matlab.obj blas\*.obj
            mex -O -largeArrayDims predict.c tron.obj linear.obj linear_model_matlab.obj blas\*.obj
            mex -O -largeArrayDims libsvmread.c
            mex -O -largeArrayDims libsvmwrite.c
        otherwise
            error('liblinear:make','Unknown platform');
    end
    cd(cwd);
end