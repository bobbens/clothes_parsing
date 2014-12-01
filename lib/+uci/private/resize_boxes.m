function boxes = resize_boxes( boxes, s )
%RESIZE_BOXES resize detections by a specified factor

for i = find(~cellfun(@isempty,boxes))
    ind = 1:(size(boxes{i},2)-2);
    boxes{i}(:,ind) = (boxes{i}(:,ind)-1) / s(i) + 1;
end

end

