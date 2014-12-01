function showskeleton(im, boxes, partcolor, components)

for i = 1:length(components),
  parent(i) = components(i).parent;
end

imagesc(im); axis image; axis off; 
if ~isempty(boxes)
  for i = 1:length(parent)
    x1(:,i) = boxes(:,1+(i-1)*4);
    y1(:,i) = boxes(:,2+(i-1)*4);
    x2(:,i) = boxes(:,3+(i-1)*4);
    y2(:,i) = boxes(:,4+(i-1)*4);
  end
  x = (x1 + x2)/2;
  y = (y1 + y2)/2;
  
  for n = 1:size(boxes,1)
    for child = 2:length(parent)
      x1 = x(n,parent(child));
      y1 = y(n,parent(child));
      x2 = x(n,child);
      y2 = y(n,child);
      line([x1 x2],[y1 y2],'color',partcolor{child},'linewidth',2);
    end
  end
end
drawnow;
