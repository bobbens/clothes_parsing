function [ cdata ] = rgbread( filename, varargin )
    %RGBREAD read an image file in rgb format
    [cdata,map] = imread( filename, varargin{:} );
    if ~isempty(map)
        cdata = ind2rgb(cdata, map);
    elseif isempty(map) && size(cdata,3)==1
        cdata = repmat(cdata,[1,1,3]);
    end
    cdata = im2uint8(cdata);
end