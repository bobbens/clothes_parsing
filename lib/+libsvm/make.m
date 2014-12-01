function make(target)
%MAKE build libsvm matlab package
%
%  libsvm.make
%  libsvm.make('all')
%  libsvm.make('clean')
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
            mex -O -largeArrayDims -c svm.cpp
            mex -O -largeArrayDims -c svm_model_matlab.c
            mex -O -largeArrayDims svmtrain.c svm.obj svm_model_matlab.obj
            mex -O -largeArrayDims svmpredict.c svm.obj svm_model_matlab.obj
            mex -O -largeArrayDims libsvmread.c
            mex -O -largeArrayDims libsvmwrite.c
        otherwise
            error('libsvm:make','Unknown platform');
    end
    cd(cwd);
end