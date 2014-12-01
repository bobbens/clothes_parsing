function [] = compile(opt)
%COMPILE build dependent mex files

setenv('MATLABDIR',matlabroot);
cwd = pwd;
root = fileparts(mfilename('fullpath'));
cd(root);
if nargin == 1 && strcmp(opt,'clean')
    cd src/gpb_src
    !make clean
    cd ../..

    cd src/buildW
    !make clean
    cd ../..
    
    cd src/savgol 
    delete(['savgol_border.',mexext]);
    cd ../..
    
    cd src/ucm
    delete(['ucm_mean_pb.',mexext]);
    cd ../..

    delete([root,'*.',mexext]);
else
    cd src/gpb_src
    !make
    % by default, mex functions use as many threads as possible
    % add NUM_THREADS=1 to limit threads
    !make matlab NUM_THREADS=1
    %!make matlab
    movefile('matlab/segmentation/mex_contour_sides.mex*',root,'f');
    movefile('matlab/segmentation/mex_nonmax_oriented.mex*',root,'f');
    movefile('matlab/segmentation/mex_pb_parts_final_selected.mex*',root,'f');
    movefile('matlab/segmentation/mex_line_inds.mex*',root,'f');
    %copyfile('matlab/segmentation/mex_contour_sides.mex*',root,'f');
    %copyfile('matlab/segmentation/mex_nonmax_oriented.mex*',root,'f');
    %copyfile('matlab/segmentation/mex_pb_parts_final_selected.mex*',root,'f');
    %copyfile('matlab/segmentation/mex_line_inds.mex*',root,'f');
    cd ../..

    cd src/buildW
    !make
    mex buildW.cpp -Iutil smatrix.cc ic.cc affinity.cc util/libutil.a
    movefile('buildW.mex*',root,'f');
    %copyfile('buildW.mex*',root,'f');
    cd ../..

    cd src/savgol 
    mex savgol_border.cpp
    movefile('savgol_border.mex*',root,'f');
    %copyfile('savgol_border.mex*',root,'f');
    cd ../..

    cd src/ucm
    mex ucm_mean_pb.cpp
    movefile('ucm_mean_pb.mex*',root,'f');
    %copyfile('ucm_mean_pb.mex*',root,'f');
    cd ../..
end
cd(cwd);

end
