function [evalNo, Score, R] = EvalStickmen(name,stickmen_set, gt_stickman, pcp_matching_threshold)
% user interface to evaluate set of stickmen together with corresponding detection bounding boxes against one ground-truth stickman
% (single-image evaluation)
%
% Input:
% stickmen_set - (contains set of stickmen to be evaluated with the ground-truth stickman) is an struct array with fields: 
%       .coor - stickman end-points coordinates (:,nparts) = [x1 y1 x2 y2]'
%       .det - the detection bounding box associated with the stickman [minx miny maxx maxy]
%
%       if .det field doesn't exist, then for association purpose a detection bounding boxes will be generated for each element in stickmen_set 
%       in the same way as for groundtruth_stickman
% gt_stickman - (contains ground-truth stickman coordiantes) is a 1x1 struct:
%        .coor - containing ground-truth stickman end-points coordinates (:,nparts) = [x1 y1 x2 y2]' 
%
% pcp_matching_threshold (optional) - defines the PCP sticks matching threshold (default 0.5) -> for definition look into README.txt 
%
% Output:
% evalNo - index of the evaluated stickman that has been matched with the ground-truth stickman
% Score -  indicates the number of correctly estimated body parts in stickmen_set(evalNo)
%
% in case when none of the stickmen_set can be associated with the gt_stickman or stickmen_set is empty then evalNo = 0, Score = nan;
%
% the evaluation procedure stops as soon one of the stickmen_set is associated to the ground-truth annotation 
% (if two det windows from stickmen_set overlap by more than PASCAL criterion then evaluation stops reporting an error)
% (at least one end-point of one segment of each stickman being evaluated must lie within the provided detection window)
%
% Eichner/2009

if nargin < 3
  pcp_matching_threshold  = 0.5;
end

Score = nan;
evalNo = 0;

R = [];

if isempty(stickmen_set) || isempty(stickmen_set(1).coor)
  return;
end

if length(gt_stickman) ~= 1 || ~isfield(gt_stickman,'coor')
  error('gt_stickman must be a 1x1 struct with field .coor')
end


N = length(stickmen_set);
classname = 'ubf';
pars.det_hwratio(2) = 0.9;
pars.iou_thresh = 0.5;


% create the .det field if does not exist
if ~isfield(stickmen_set(1),'det')
 stickmen_set(1).det = []; % setting for one elements adds this field for whole array
end
for i=1:N
  if isempty(stickmen_set(i).det)
    % if detection was not provided then derive from stickman
    if strcmp(name,'PARSE')
      stickmen_set(i).det = PARSE_detBBFromStickman(stickmen_set(i).coor);
    else
      stickmen_set(i).det = detBBFromStickman(stickmen_set(i).coor,classname,pars);
    end
  end
end

% no two of the provided detections may be overlaping by more than the PASCAL criterion
% overlaping = zeros(N);
% for i=1:N
%   for j=(i+1):N
%     overlaping(i,j) = all_pairs_bb_iou(stickmen_set(i).det', stickmen_set(j).det');
%   end
% end
% if any(overlaping(:) > pars.iou_thresh)
%   error('no two detections may overlap with each other by more than the PASCAL criterion !!!!');
% end

% hallucinate detection window from the gt stickman
if strcmp(name,'PARSE')
  gt_bb = PARSE_detBBFromStickman(gt_stickman.coor);
else
  gt_bb = detBBFromStickman(gt_stickman.coor,classname,pars);
end

for i=1:N
  if all_pairs_bb_iou(gt_bb',  stickmen_set(i).det') > pars.iou_thresh 
    % evaluate only when eval_bb can be associated with gt_bb iou > 0.5
    % stop evaluation after first bb match
    evalNo = i; % so evalNo ~= 0 and stickman is taken into account when calculating PCP
    % check whether a stickman truly lie in its detection window 
    %(at least single end-point of one segment must lie inside the corresponding bb)
    if ~any(boxTest(stickmen_set(i).det, [stickmen_set(i).coor(1:2,:) stickmen_set(i).coor(3:4,:)]))
      disp('WARNING: evaluated stickmen is not inside its detection bounding box!!!! -> score 0');
      Score = 0;
    else
      [Score R] = DirectEvalStickman(stickmen_set(i).coor, gt_stickman.coor,pcp_matching_threshold);
    end
    return
  end
end

end

function bb = detBBFromStickman(stick,classname,params)
% function minBB(stick,classname)
% stick in format rows: [x1; y1; x2; y2] cols: nsticks
% returns bb in format [minx miny maxx maxy]
% so far just for ubf
% requires params.det_hwratio(classid)
stick = double(stick);
classid = class_name2id(classname);
switch classname
  case 'ubf'
    %minx = min([stick(1,[ 1 2 3 6]) stick(3,[1 2 3 6])]); % min x of all parts besides lower arms
    %maxx = max([stick(1,[2 3]) stick(3,[2 3])]); % max x of all parts besides lower arms
    torso_center = [(stick(1,1)+stick(3,1))/2 (stick(2,1)+stick(4,1))/2]; %center [x y]
    miny = min([stick(2,[1 6]) stick(4,[1 6]) torso_center(2)]); % min y of torso and head
    maxy = max([stick(2,6) stick(4,6) torso_center(2) ]); %max y of head and torso center
    diffx = abs(miny-maxy)/params.det_hwratio(classid);
    head_center = [(stick(1,6)+stick(3,6))/2 (stick(2,6)+stick(4,6))/2];
    minx = (head_center(1)+torso_center(1))/2 - diffx/2;
    maxx = (head_center(1)+torso_center(1))/2 + diffx/2;
  case 'ubp'
    error('not written yet');
  case 'full'
    error('not written yet');
end
bb = [minx miny maxx maxy];
end

function id = class_name2id(name)
% id = class_name2id(name)
%
% object class name-to-id convertor
%
% names are case-insensitive
%
% See also class_id2name
%

switch lower(name)
    case 'ubp'
        id = 1;
    case 'ubf'
        id = 2;
    case 'full'
        id = 3;
    case 'ubfreg' %if it is from regressed ubf from face detection treat the same as normal ubf 
        id = 2;
    otherwise
        error([mfilename ': unknown class name ' name]);
end
end

function bb = PARSE_detBBFromStickman(stick)
stick = double(stick);

torso_center = [(stick(1,1)+stick(3,1))/2 (stick(2,1)+stick(4,1))/2]; %center [x y]
miny = min([stick(2,[1 10]) stick(4,[1 10]) torso_center(2)]); % min y of torso and head
maxy = max([stick(2,10) stick(4,10) torso_center(2)]); % max y of head and torso center
diffx = abs(miny-maxy)/0.9;
head_center = [(stick(1,10)+stick(3,10))/2 (stick(2,10)+stick(4,10))/2];
minx = (head_center(1)+torso_center(1))/2 - diffx/2;
maxx = (head_center(1)+torso_center(1))/2 + diffx/2;
bb = [minx miny maxx maxy];
end

function test = boxTest(box,points)
% boxTest(box,points)
% test weather point lies in the box
% Input:
%  point(:,n) = [x y]'
%  box = [minx miny maxx maxy]
% Output:
%   test = bool(1,n) - true/false value for each point passed to the routine
%
points = points';
test = (points(:,1) > box(1) & points(:,1) < box(3)) | (points(:,2) > box(2) & points(:,2) < box(4));
test = test';
end


function [Score R] = DirectEvalStickman(estimated_stickman,groundtruth_stickman, pcp_matching_threshold)
% user interface to evaluate single stickman against a single ground-truth stickman
% estimated_stickman and groundtruth_stickman are arrays of stickman end-points coordinates (:,nparts) = [x1 y1 x2 y2]'
% pcp_matching_threshold (optional) - defines the PCP sticks matching threshold (default 0.5) -> for definition look into README.txt 
%
% Output:
% Score - number of correctly estimated body parts 
% R(p) - indicator of parts that have been correctly estimated
%
% Eichner/2009

if nargin < 3
  pcp_matching_threshold = 0.5;
end

% evaluation criterion used in Ferrari (cvpr08, cvpr09) and Eichner (bmvc09) papers
pars.max_parall_dist = pcp_matching_threshold;
pars.max_perp_dist = pcp_matching_threshold;

R = DirectEvalSegms(double(estimated_stickman),double(groundtruth_stickman),pars);
Score = sum(R);
end

function R = DirectEvalSegms(S1, S2, pars)

% Compare segments S1 to S2.
%
% -> S(:,six) = [x1 y1 x2 y2]'
%
% if pars.crit == 'col'
% Two segments are equivalent if
% within the thresholds in pars.max_*
%
% if pars.crit == 'endpts'
% Two segmentas are equivalent if
% their endpts are within pars.max_{parall,perp}_dist
%
% Output:
% - R(six) = true iff S1(:,six) equivalent to S2(:,six)
%

S1_col = XY2COL(S1);   
S2_col = XY2COL(S2);

Np = size(S2,2);

%norm factor for classification
norm_fact = S2_col(4,:);

for p = 1:Np
  %
  % current segment1 and segment2
  cs1 = S1(:,p);
  cs2 = S2(:,p);
  %
  % match endpts so that order is same in cs1 and cs2
  temp = vgg_nearest_neighbour_dist(cs1(1:2), reshape(cs2,2,2));
  [trash flip] = min(temp);
  if flip == 2
    cs1 = cs1([3 4 1 2]);
  end
  %
  % distances
  dxy = cs1 - cs2;
  ndxy = zeros(4,1);
  % rotate onto cs2 direction
  theta = -S2_col(3,p);                          % must rotate back
  Rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
  ndxy(1:2) = Rot * dxy(1:2);                    % now dx -> distance parallel to S2 dir, and dy -> dist perpendicular to it
  ndxy(3:4) = Rot * dxy(3:4);      
  %
  % normalize
  ndxy = ndxy / norm_fact(p);
  %
  % classify
  ndxy = abs(ndxy);
  avg_ndxy = [mean(ndxy([1 3])) mean(ndxy([2 4]))];
  R(p) = (avg_ndxy(1) < pars.max_parall_dist) && (avg_ndxy(2) < pars.max_perp_dist);
  %if p == 4
  %    keyboard;
  %end
end

end

function Scol = XY2COL(S)

% converts segments S(:,ix) = [x1 y1 x2 y2]'
% to Scol(:,ix) = [ctr_x ctr_y orient length]'
%
% with orient in [0,pi]
%

ctrs = [mean(S([1 3], :)); mean(S([2 4], :))];
dx = S(3,:)-S(1,:);
dy = S(4,:)-S(2,:);
lens = sqrt(dx.^2+dy.^2);
orients = atan(dy./dx);  % in [-pi/2,pi/2]
neg = orients < 0;
orients(neg) = orients(neg) + pi;

Scol = [ctrs; orients; lens];
end