TARGETDIR=../../

cd gpb_src
make clean
make
make matlab
cp -f matlab/segmentation/mex_contour_sides.mex* $TARGETDIR
cp -f matlab/segmentation/mex_nonmax_oriented.mex* $TARGETDIR
cp -f matlab/segmentation/mex_pb_parts_final_selected.mex* $TARGETDIR
cp -f matlab/segmentation/mex_line_inds.mex* $TARGETDIR
cd ..

cd buildW
make clean
make
cp -f buildW.mex* $TARGETDIR
cd ..

cd savgol 
matlab -nodisplay -nojvm -r "mex savgol_border.cpp; exit"
cp -f savgol_border.mex* $TARGETDIR
cd ..

cd ucm
matlab -nodisplay -nojvm -r "mex ucm_mean_pb.cpp; exit"
cp -f ucm_mean_pb.mex* $TARGETDIR
cd ..
