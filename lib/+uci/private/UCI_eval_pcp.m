function [detRate PCP R] = UCI_eval_pcp(name,boxes,test)
% Stick parts order:
% torso, ul_leg, ur_leg, ll_leg, lr_leg, ul_arm, ur_arm, ll_arm, lr_arm, head

name = 'PARSE';

% Convert boxes to points
points = cell(size(boxes));
for i = 1:numel(boxes)
  if ~isempty(boxes{i})
    for p = 1:floor(size(boxes{i},2)/4),
      cx1 = boxes{i}(:,1+(p-1)*4);
      cy1 = boxes{i}(:,2+(p-1)*4);
      cx2 = boxes{i}(:,3+(p-1)*4);
      cy2 = boxes{i}(:,4+(p-1)*4);

      points{i}(:,1+(p-1)*2) = (cx1+cx2)/2;
      points{i}(:,2+(p-1)*2) = (cy1+cy2)/2;
    end
  end
end

% -------------------
% generate testing stick
I = [1   1   2   2   3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20];
J = [3   15  10  22  10 12 22 24 12 14 24 26 3  5  15 17 5  7  17 19 1  2];
S = [1/2 1/2 1/2 1/2 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1];
A = full(sparse(I,J,S,20,26));

for n = 1:length(test)
  pts = points{n};
  if isempty(pts)
    predstick(n).stickmen = [];
    continue;
  end    
  for i = 1:size(pts,1)
    predstick(n).stickmen(i).coor = reshape((A*reshape(pts(i,:),2,size(pts,2)/2)')',4,10);
  end
end

% -------------------
% create groundtruth stick
% because PARSE dataset do not have grountruth stick labels, we create the
% groundtruth stick labels by ourselves using groundtruth keypoints

% UCI_TO_PARSE = full(sparse(...
%     [14;13;9;8;7;3;2;1;10;11;12;4;5;6],...
%     [1;2;3;5;7;10;12;14;15;17;19;22;24;26],...
%     1, 14, 26));
% 
% I = [1   1   2   2   3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
% J = [9   10  3   4   3 2 4 5 2 1 5 6  9  8  10 11 8  7  11 12 14 13];
% S = [1/2 1/2 1/2 1/2 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1  1  1  1];
% A = full(sparse(I,J,S,20,14)) * UCI_TO_PARSE;

A = full(sparse(...
    [19,20,1,11,12,15,16,2,3,4,7,8,1,13,14,17,18,2,5,6,9,10],...
    [1,2,3,3,5,5,7,10,10,12,12,14,15,15,17,17,19,22,22,24,24,26],...
    [1,1,0.5,1,1,1,1,0.5,1,1,1,1,0.5,1,1,1,1,0.5,1,1,1,1], 20, 26));

for n = 1:length(test)
  gtstick(n).stickmen.coor = reshape((A*test(n).point)',4,10);
end

% the PCP evaluation function originally comes from BUFFY dataset, we keep using that for performance evaluation
[detRate PCP R] = eval_pcp(name,predstick,gtstick);
