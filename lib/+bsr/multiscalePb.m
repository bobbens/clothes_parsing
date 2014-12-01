function [mPb_nmax, mPb_nmax_rsz, bg1, bg2, bg3, cga1, cga2, cga3, cgb1, cgb2, cgb3, tg1, tg2, tg3, textons] = multiscalePb(im, rsz)
%function [mPb_nmax, mPb_nmax_rsz, bg1, bg2, bg3, cga1, cga2, cga3, cgb1, cgb2, cgb3, tg1, tg2, tg3, textons] = multiscalePb(im, rsz)
%
% description:
% compute local contour cues of an image.
%
% gradients by Michael Maire <mmaire@eecs.berkeley.edu>
%
% Pablo Arbelaez <arbelaez@eecs.berkeley.edu>
% December 2010
if nargin<2, rsz = 1.0; end


% default feature weights
if size(im,3) == 3,
    weights = [0.0146    0.0145    0.0163    0.0210    0.0243    0.0287    0.0166    0.0185    0.0204    0.0101    0.0111    0.0141];
else
    im(:,:,2)=im(:,:,1);im(:,:,3)=im(:,:,1);
    weights = [0.0245    0.0220    0.0233         0         0         0         0         0         0    0.0208    0.0210    0.0229];
end

% get gradients
%tic;
[bg1, bg2, bg3, cga1, cga2, cga3, cgb1, cgb2, cgb3, tg1, tg2, tg3, textons] = det_mPb(im);
%fprintf('Local cues: %g\n', toc);

% smooth cues
gtheta = [1.5708    1.1781    0.7854    0.3927   0    2.7489    2.3562    1.9635];
% tic;
filters = make_filters([3 5 10 20], gtheta);
for o = 1 : size(tg1, 3),
    bg1(:,:,o) = fitparab(bg1(:,:,o),3,3/4,gtheta(o),filters{1,o});
    bg2(:,:,o) = fitparab(bg2(:,:,o),5,5/4,gtheta(o),filters{2,o});
    bg3(:,:,o) = fitparab(bg3(:,:,o),10,10/4,gtheta(o),filters{3,o});

    cga1(:,:,o) = fitparab(cga1(:,:,o),5,5/4,gtheta(o),filters{2,o});
    cga2(:,:,o) = fitparab(cga2(:,:,o),10,10/4,gtheta(o),filters{3,o});
    cga3(:,:,o) = fitparab(cga3(:,:,o),20,20/4,gtheta(o),filters{4,o});

    cgb1(:,:,o) = fitparab(cgb1(:,:,o),5,5/4,gtheta(o),filters{2,o});
    cgb2(:,:,o) = fitparab(cgb2(:,:,o),10,10/4,gtheta(o),filters{3,o});
    cgb3(:,:,o) = fitparab(cgb3(:,:,o),20,20/4,gtheta(o),filters{4,o});

    tg1(:,:,o) = fitparab(tg1(:,:,o),5,5/4,gtheta(o),filters{2,o});
    tg2(:,:,o) = fitparab(tg2(:,:,o),10,10/4,gtheta(o),filters{3,o});
    tg3(:,:,o) = fitparab(tg3(:,:,o),20,20/4,gtheta(o),filters{4,o});

end
% fprintf('Cues smoothing: %g\n', toc);


% compute mPb at full scale
mPb_all = zeros(size(tg1));
for o = 1 : size(mPb_all, 3),
    l1 = weights(1)*bg1(:, :, o);
    l2 = weights(2)*bg2(:, :, o);
    l3 = weights(3)*bg3(:, :, o);

    a1 = weights(4)*cga1(:, :, o);
    a2 = weights(5)*cga2(:, :, o);
    a3 = weights(6)*cga3(:, :, o);

    b1 = weights(7)*cgb1(:, :, o);
    b2 = weights(8)*cgb2(:, :, o);
    b3 = weights(9)*cgb3(:, :, o);

    t1 = weights(10)*tg1(:, :, o);
    t2 = weights(11)*tg2(:, :, o);
    t3 = weights(12)*tg3(:, :, o);

    mPb_all(:, :, o) = l1 + a1 + b1 + t1 + l2 + a2 + b2 + t2 + l3 + a3 + b3 + t3;

end

% non-maximum suppression
mPb_nmax = nonmax_channels(mPb_all);
mPb_nmax = max(0, min(1, 1.2*mPb_nmax));


% compute mPb_nmax resized if necessary
if rsz < 1,
    mPb_all = imresize(tg1, rsz);
    mPb_all(:) = 0;

    for o = 1 : size(mPb_all, 3),
        l1 = weights(1)*bg1(:, :, o);
        l2 = weights(2)*bg2(:, :, o);
        l3 = weights(3)*bg3(:, :, o);

        a1 = weights(4)*cga1(:, :, o);
        a2 = weights(5)*cga2(:, :, o);
        a3 = weights(6)*cga3(:, :, o);

        b1 = weights(7)*cgb1(:, :, o);
        b2 = weights(8)*cgb2(:, :, o);
        b3 = weights(9)*cgb3(:, :, o);

        t1 = weights(10)*tg1(:, :, o);
        t2 = weights(11)*tg2(:, :, o);
        t3 = weights(12)*tg3(:, :, o);

        mPb_all(:, :, o) = imresize(l1 + a1 + b1 + t1 + l2 + a2 + b2 + t2 + l3 + a3 + b3 + t3, rsz);

    end

    mPb_nmax_rsz = nonmax_channels(mPb_all);
    mPb_nmax_rsz = max(0, min(1, 1.2*mPb_nmax_rsz));
else
    mPb_nmax_rsz = mPb_nmax;
end
end

%%
function filters = make_filters(radii, gtheta)

d = 2; 

filters = cell(numel(radii), numel(gtheta));
for r = 1:numel(radii),
    for t = 1:numel(gtheta),
        
        ra = radii(r);
        rb = ra / 4;
        theta = gtheta(t);
        
        ra = max(1.5, ra);
        rb = max(1.5, rb);
        ira2 = 1 / ra^2;
        irb2 = 1 / rb^2;
        wr = floor(max(ra, rb));
        wd = 2*wr+1;
        sint = sin(theta);
        cost = cos(theta);
        
        % 1. compute linear filters for coefficients
        % (a) compute inverse of least-squares problem matrix
        filt = zeros(wd,wd,d+1);
        xx = zeros(2*d+1,1);
        for u = -wr:wr,
            for v = -wr:wr,
                ai = -u*sint + v*cost; % distance along major axis
                bi = u*cost + v*sint; % distance along minor axis
                if ai*ai*ira2 + bi*bi*irb2 > 1, continue; end % outside support
                xx = xx + cumprod([1;ai+zeros(2*d,1)]);
            end
        end
        A = zeros(d+1,d+1);
        for i = 1:d+1,
            A(:,i) = xx(i:i+d);
        end
        
        % (b) solve least-squares problem for delta function at each pixel
        for u = -wr:wr,
            for v = -wr:wr,
                ai = -u*sint + v*cost; % distance along major axis
                bi = u*cost + v*sint; % distance along minor axis
                if (ai*ai*ira2 + bi*bi*irb2) > 1, continue; end % outside support
                yy = cumprod([1;ai+zeros(d,1)]);
                filt(v+wr+1,u+wr+1,:) = A\yy;
            end
        end
        
        filters{r,t}=filt;
    end
end
end

function [bg1, bg2, bg3, cga1, cga2, cga3, cgb1, cgb2, cgb3, tg1, tg2, tg3, textons] = det_mPb(im)
% compute image gradients. Implementation by Michael Maire.

% compute pb parts
[ ...
    textons, ...
    bg_r3, bg_r5,  bg_r10,  cga_r5, cga_r10, cga_r20, cgb_r5, cgb_r10, cgb_r20, tg_r5,  tg_r10,  tg_r20...
    ] = bsr.mex_pb_parts_final_selected(im(:,:,1),im(:,:,2),im(:,:,3));

[sx sy sz] = size(im);
temp = zeros([sx sy 8]);

for r = [3 5 10]
    for ori = 1:8
        eval(['temp(:,:,ori) = bg_r' num2str(r) '{' num2str(ori) '};']);
    end
    eval(['bg_r' num2str(r) ' = temp;']);
end
bg1 = bg_r3; bg2 = bg_r5;  bg3 = bg_r10; 

for r = [5 10 20]
    for ori = 1:8
        eval(['temp(:,:,ori) = cga_r' num2str(r) '{' num2str(ori) '};']);
    end
    eval(['cga_r' num2str(r) ' = temp;']);
end
cga1 = cga_r5; cga2 = cga_r10;  cga3 = cga_r20; 

for r = [5 10 20]
    for ori = 1:8
        eval(['temp(:,:,ori) = cgb_r' num2str(r) '{' num2str(ori) '};']);
    end
    eval(['cgb_r' num2str(r) ' = temp;']);
end
cgb1 = cgb_r5; cgb2 = cgb_r10;  cgb3 = cgb_r20; 

for r = [5 10 20]
    for ori = 1:8
        eval(['temp(:,:,ori) = tg_r' num2str(r) '{' num2str(ori) '};']);
    end
    eval(['tg_r' num2str(r) ' = temp;']);
end
tg1 = tg_r5; tg2 = tg_r10;  tg3 = tg_r20; 
end

%%
function a = fitparab(z,ra,rb,theta,filt)
% function [a,b,c] = fitparab(z,ra,rb,theta)
%
% Fit cylindrical parabolas to elliptical patches of z at each
% pixel.  
%
% INPUT
%	z	Values to fit.
%	ra,rb	Radius of elliptical neighborhood, ra=major axis.
%	theta	Orientation of fit (i.e. of minor axis).
%
% OUTPUT
%	a,b,c	Coefficients of fit: a + bx + cx^2
%


% compute the interior quickly with convolutions
a = conv2(z,filt(:,:,1),'same');
%fix border with mex file
a = bsr.savgol_border(a, z, ra, rb, theta);
end

%%
function nmax = nonmax_channels(pb, nonmax_ori_tol)
% given NxMxnum_ori oriented channels, compute oriented nonmax suppression
if (nargin < 2), nonmax_ori_tol = pi/8; end
n_ori = size(pb,3);
oris = (0:(n_ori-1))/n_ori * pi;
[y,i] = max(pb,[],3);
i = oris(i);
y(y<0)=0;
nmax = nonmax_oriented(y,i, nonmax_ori_tol);
end

%%
function nmax = nonmax_oriented(pb, ori, nonmax_ori_tol)
% Oriented non-max suppression (2D).
%
% Perform non-max suppression orthogonal to the specified orientation on
% the given 2D matrix using linear interpolation in a 3x3 neighborhood.
%
% A local maximum must be greater than the interpolated values of its 
% adjacent elements along the direction orthogonal to this orientation.
%
% If an orientation is specified per element, then the elements themselves
% may optionally be treated as oriented vectors by specifying a value less 
% than pi/2 for the orientation tolerance.  In this case, neighboring 
% vectors are projected along a line in the local orientation direction and
% the lengths of the projections are used in determining local maxima.
% When projecting, the orientation tolerance is subtracted from the true
% angle between the vector and the line (with a result less than zero
% causing the length of the projection to be the length of the vector).
%
% Non-max elements are assigned a value of zero.
%
% NOTE: The original matrix must be nonnegative.
if (nargin < 3), nonmax_ori_tol = pi/8; end
nmax = bsr.mex_nonmax_oriented(pb, ori, nonmax_ori_tol);
end
