function [boxes] = testmodel(name,model,X,suffix,verbose)

%globals;
if nargin < 5, verbose = true; end

% try
%   load([cachedir name '_boxes_' suffix]);
% catch
  boxes = cell(1,length(X));
%   points = cell(1,length(X));
  for i = 1:length(X)
    if verbose
        fprintf([name ': testing: %d/%d\n'],i,length(X));
    end
    %im = imread(X(i).im);
    %BEGIN HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    box = detect(X(i),model,model.thresh);
    %END HACK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %box = detect(im,model,model.thresh);
    if ~isempty(box)
      boxes{i} = nms(box,0.5);
%       for p = 1:floor(size(boxes{i},2)/4),
%         cx1 = boxes{i}(:,1+(p-1)*4);
%         cy1 = boxes{i}(:,2+(p-1)*4);
%         cx2 = boxes{i}(:,3+(p-1)*4);
%         cy2 = boxes{i}(:,4+(p-1)*4);
%         
%         points{i}(:,1+(p-1)*2) = (cx1+cx2)/2;
%         points{i}(:,2+(p-1)*2) = (cy1+cy2)/2;
%       end
    else
      boxes{i} = [];
%       points{i} = [];
    end
  end

  if nargin < 3
    suffix = [];
  end
%   save([cachedir name '_boxes_' suffix], 'boxes','points','model');
% end
