% load image
im = double(imread('../../test/io/formats/image/jpeg_test/images/test_color.jpg'))/255;

% compute pb
[bg, cg, tg, pb, pb_nmax] = mex_pb(im(:,:,1),im(:,:,2),im(:,:,3));

% threshold pb
pb_nmax_thresh = pb_nmax .* (pb_nmax > 0.1);

% extract contours
contours = contour(pb_nmax_thresh);

% display contours
disp_contours(contours, im);
