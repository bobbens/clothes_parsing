function theta = std_angle(theta, min_angle)
% STD_ANGLE    Standardize a matrix of angle values to a specified range.
%
%  THETA = STD_ANGLE(theta) adjusts angle values to the [-pi,pi) range.
%
%  THETA = STD_ANGLE(theta, min_angle) adjusts angle values to be in the
%  range [min_angle, min_angle + 2*pi).

% check arguments
error(nargchk(1,2,nargin));
if (nargout > 1), error('Too many output arguments.'); end

% set default arguments
if (nargin < 2), min_angle = -pi; end
max_angle = min_angle + 2*pi;

% adjust angle values that are too small
small_inds = find(theta < min_angle);
theta(small_inds) = theta(small_inds) + ceil((min_angle - theta(small_inds))/(2*pi))*(2*pi);

% adjust angle values that are too big
big_inds = find(theta > max_angle);
theta(big_inds) = theta(big_inds) - ceil((theta(big_inds) - max_angle)/(2*pi))*(2*pi);

% adjust angle values on edge of range
theta(find(theta == max_angle)) = min_angle;
