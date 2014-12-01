function plot_points(x, y, ori, scale, N, varargin)
% PLOT_POINTS  Plot oriented polygons centered at specified points.
%
%  PLOT_POINTS(x, y, ori, scale, N) plots N-sided polygons at the given 
%  locations, scales, and orientations.  The scale specifies the radius
%  of the polygon's circumcircle.  The scale, ori, and N arguments may 
%  be scalars or vectors (specifying global or per-polygon parameters).
%
%  PLOT_POINTS(x, y, ...) where some or all additional parameters are 
%  omitted uses default values for those missing parameters (defaults are
%  scale=1, ori=pi/2, N=3).
%
%  PLOT_POINTS(y) uses default parameters and sets x = 1:length(y).
%
%  PLOT_POINTS(x, y, ori, scale, N, ...) passes additional arguments to 
%  Matlab's plot commmand, which is called to draw the line segments for
%  the polygons.

% check arguments
error(nargchk(1,Inf,nargin));
if (nargout > 0), error('Too many output arguments.'); end

% set default arguments
n_points = length(x);
if (nargin < 2), y = x; x = 1:n_points; end
if (nargin < 3), ori = (pi/2).*ones([1 n_points]); end
if (nargin < 4), scale = ones([1 n_points]); end
if (nargin < 5), N = 8.*ones([1 n_points]); end

% check argument lengths
if (isscalar(ori)), ori = repmat(ori,[1 n_points]); end
if (isscalar(scale)), scale = repmat(scale,[1 n_points]); end
if (isscalar(N)), N = repmat(N, [1 n_points]); end

% save hold state
hold_save = ishold;

% plot each point
for n = 1:n_points
   % corners of polygon
   px = [cos((0:N(n))/N(n)*(2*pi)) 0];
   py = [sin((0:N(n))/N(n)*(2*pi)) 0];
   % rotate polygon by descriptor ori
   % and scale by descriptor size
   rot_mx = scale(n).*[ ...   
       cos(ori(n)) -sin(ori(n)); ...
       sin(ori(n))  cos(ori(n)) ...
   ];
   p = rot_mx*[px; py];
   px = p(1,:) + x(n);
   py = p(2,:) + y(n);
   % plot the polygon
   plot(px, py, varargin{:});
   % hold after the first point is plotted
   hold on;
end

% restore hold state
if (hold_save == 0)
    hold off;
end
