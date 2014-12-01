function [sPb] = spectralPb(mPb, orig_sz, nvec)
% function [sPb] = spectralPb(mPb, orig_sz, outFile, nvec)
%
% description:
%   global contour cue from local mPb.
%
% computes Intervening Contour with BSE code by Charless Fowlkes:
%
%http://www.cs.berkeley.edu/~fowlkes/BSE/
%
% Pablo Arbelaez <arbelaez@eecs.berkeley.edu>
% December 2010

if nargin<3, nvec = 17; end
if nargin<2, orig_sz = size(mPb); end

[tx, ty] = size(mPb);

l{1} = zeros(size(mPb, 1) + 1, size(mPb, 2));
l{1}(2:end, :) = mPb;
l{2} = zeros(size(mPb, 1), size(mPb, 2) + 1);
l{2}(:, 2:end) = mPb;

% build the pairwise affinity matrix
[val,I,J] = bsr.buildW(l{1},l{2});
W = sparse(val,I,J);

[wx, wy] = size(W);
x = 1 : wx;
S = full(sum(W, 1));
D = sparse(x, x, S, wx, wy);
clear S x;

opts.issym=1;
opts.isreal = 1;
%opts.disp=2;

% Try eigenvalue decomposition
i = nvec;
sigma = 0;
while i >= 1
    try
        [EigVect, EVal] = eigs(D - W, D, i, sigma, opts);
        break;
    catch e
        disp(e.getReport);
        i = floor(i / 2);
        if sigma==0
            sigma = eps;
	else
            sigma = 2*sigma;
        end
        opts.disp=2;
    end
end
nvec = i;
if nvec < 1 || ~exist('EigVect','var')
    error('bsr:spectralPb:badInput','The input is singular');
end
clear D W opts i;

EigVal = diag(EVal);
clear Eval;

EigVal(1:end) = EigVal(end:-1:1);
EigVect(:, 1:end) = EigVect(:, end:-1:1);

txo=orig_sz(1); tyo=orig_sz(2); 
vect = zeros(txo, tyo, nvec);
for v = 2 : nvec,
    vect(:, :, v) = imresize(reshape(EigVect(:, v), [ty tx])',[txo,tyo]);
end
clear EigVect;

% spectral Pb
for v=2:nvec,
    vect(:,:,v)=(vect(:,:,v)-min(min(vect(:,:,v))))/(max(max(vect(:,:,v)))-min(min(vect(:,:,v))));
end

% OE parameters
hil = 0;
deriv = 1;
support = 3;
sigma = 1;
norient = 8;
dtheta = pi/norient;
ch_per = [4 3 2 1 8 7 6 5];

sPb = zeros(txo, tyo, norient);
for v = 1 : nvec
    if EigVal(v) > 0,
        vec = vect(:,:,v)/sqrt(EigVal(v));
        for o = 1 : norient,
            theta = dtheta*o;
            f = oeFilter(sigma, support, theta, deriv, hil);
            sPb(:,:,ch_per(o)) = sPb(:,:,ch_per(o)) + abs(applyFilter(f, vec));
        end
    end
end

end

%%
function [f] = oeFilter(sigma,support,theta,deriv,hil,vis)
% function [f] = oeFilter(sigma,support,theta,deriv,hil,vis)
%
% Compute unit L1-norm 2D filter.
% The filter is a Gaussian in the x direction.
% The filter is a Gaussian derivative with optional Hilbert
% transform in the y direction.
% The filter is zero-meaned if deriv>0.
%
% INPUTS
%	sigma		Scalar, or 2-element vector of [sigmaX sigmaY].
%	[support]	Make filter +/- this many sigma.
%	[theta]		Orientation of x axis, in radians.
%	[deriv]		Degree of y derivative, one of {0,1,2}.
%	[hil]		Do Hilbert transform in y direction?
%	[vis]		Visualization for debugging?
%
% OUTPUTS
%	f	Square filter.
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% March 2003

nargchk(1,6,nargin);
if nargin<2, support=3; end
if nargin<3, theta=0; end
if nargin<4, deriv=0; end
if nargin<5, hil=0; end
if nargin<6, vis=0; end

if numel(sigma)==1,
  sigma = [sigma sigma];
end
if deriv<0 | deriv>2,
  error('deriv must be in [0,2]');
end

% Calculate filter size; make sure it's odd.
hsz = max(ceil(support*sigma));
sz = 2*hsz + 1;

% Sampling limits.
maxsamples = 1000; % Max samples in each dimension.
maxrate = 10; % Maximum sampling rate.
frate = 10; % Over-sampling rate for function evaluation.

% Cacluate sampling rate and number of samples.
rate = min(maxrate,max(1,floor(maxsamples/sz)));
samples = sz*rate;

% The 2D samping grid.
r = floor(sz/2) + 0.5 * (1 - 1/rate);
dom = linspace(-r,r,samples);
[sx,sy] = meshgrid(dom,dom);

% Bin membership for 2D grid points.
mx = round(sx);
my = round(sy);
membership = (mx+hsz+1) + (my+hsz)*sz;

% Rotate the 2D sampling grid by theta.
su = sx*sin(theta) + sy*cos(theta);
sv = sx*cos(theta) - sy*sin(theta);

if vis,
  figure(1); clf; hold on;
  plot(sx,sy,'.');
  plot(mx,my,'o');
  %plot([sx(:) mx(:)]',[sy(:) my(:)]','k-');
  plot(su,sv,'x');
  axis equal;
  ginput(1);
end

% Evaluate the function separably on a finer grid.
R = r*sqrt(2)*1.01; % radius of domain, enlarged by >sqrt(2)
fsamples = ceil(R*rate*frate); % number of samples
fsamples = fsamples + mod(fsamples+1,2); % must be odd
fdom = linspace(-R,R,fsamples); % domain for function evaluation
gap = 2*R/(fsamples-1); % distance between samples

% The function is a Gaussian in the x direction...
fx = exp(-fdom.^2/(2*sigma(1)^2));
% .. and a Gaussian derivative in the y direction...
fy = exp(-fdom.^2/(2*sigma(2)^2));
switch deriv,
 case 1,
  fy = fy .* (-fdom/(sigma(2)^2));
 case 2,
  fy = fy .* (fdom.^2/(sigma(2)^2) - 1);
end
% ...with an optional Hilbert transform.
if hil,
  fy = imag(hilbert(fy));
end

% Evaluate the function with NN interpolation.
xi = round(su/gap) + floor(fsamples/2) + 1;
yi = round(sv/gap) + floor(fsamples/2) + 1;
f = fx(xi) .* fy(yi);

% Accumulate the samples into each bin.
f = isum(f,membership,sz*sz);
f = reshape(f,sz,sz);

% zero mean
if deriv>0,
  f = f - mean(f(:));
end

% unit L1-norm
sumf = sum(abs(f(:)));
if sumf>0,
  f = f / sumf;
end
end

%%
function acc = isum(x,idx,nbins)
% function acc = isum(x,idx,nbins)
%
% Indexed sum reduction, where acc(i) contains the sum of
% x(find(idx==i)).
%
% The mex version is 300x faster in R12, and 4x faster in R13.  As far
% as I can tell, there is no way to do this efficiently in matlab R12.
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% March 2003

acc = zeros(nbins,1);
for i = 1:numel(x),
  if idx(i)<1, continue; end
  if idx(i)>nbins, continue; end
  acc(idx(i)) = acc(idx(i)) + x(i);
end
end

%%
function [fim] = applyFilter(f,im)
% function [fim] = applyFilter(f,im)
%
% Apply a filter to an image with reflected boundary conditions.
%
% See also fbCreate, fbRun.
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% March 2003

fim = fbRun({f},im);
fim = fim{1};
end

%%
function [fim] = fbRun(fb,im)
% function [fim] = fbRun(fb,im)
%
% Run a filterbank on an image with reflected boundary conditions.
%
% See also fbCreate,padReflect.
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% March 2003

% find the max filter size
maxsz = max(size(fb{1}));
for i = 1:numel(fb),
  maxsz = max(maxsz,max(size(fb{i})));
end

% pad the image 
r = floor(maxsz/2);
impad = padReflect(im,r);

% run the filterbank on the padded image, and crop the result back
% to the original image size
fim = cell(size(fb));
for i = 1:numel(fb),
  if size(fb{i},1)<50,
    fim{i} = conv2(impad,fb{i},'same');
  else
    fim{i} = fftconv2(impad,fb{i});
  end
  fim{i} = fim{i}(r+1:end-r,r+1:end-r);
end
end

%%
function [impad] = padReflect(im,r)
% function [impad] = padReflect(im,r)
%
% Pad an image with a border of size r, and reflect the image into
% the border.
%
% David R. Martin <dmartin@eecs.berkeley.edu>
% March 2003

impad = zeros(size(im)+2*r);
impad(r+1:end-r,r+1:end-r) = im; % middle
impad(1:r,r+1:end-r) = flipud(im(1:r,:)); % top
impad(end-r+1:end,r+1:end-r) = flipud(im(end-r+1:end,:)); % bottom
impad(r+1:end-r,1:r) = fliplr(im(:,1:r)); % left
impad(r+1:end-r,end-r+1:end) = fliplr(im(:,end-r+1:end)); % right
impad(1:r,1:r) = rot90(im(1:r,1:r),2); % top-left
impad(1:r,end-r+1:end) = rot90(im(1:r,end-r+1:end),2); % top-right
impad(end-r+1:end,1:r) = rot90(im(end-r+1:end,1:r),2); % bottom-left
impad(end-r+1:end,end-r+1:end) = rot90(im(end-r+1:end,end-r+1:end),2); % bottom-right
end
