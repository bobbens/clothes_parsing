classdef BasePhoto < handle
    %BASEPHOTO abstract class for data instances in clothing parser
    
    properties (Abstract)
        id                % instance id
        width             % width of an image
        height            % height of an image
        image_path        % path to an image file (or in-memory data)
        segmentation      % segmentation object
        pose              % pose object
    end
    
    properties (Dependent)
        image             % image loader
    end
    
    methods
        function im = get.image(this)
            %GET.IMAGE load image from a file
            p = this.image_path;
            if ischar(this.image_path) && exist(p,'file')
                im = imread(p);
            elseif isnumeric(this.image_path)
                im = this.image_path;
            else
                im = [];
            end
            if (~isempty(this.height) && this.height ~= size(im,1)) ||...
                (~isempty(this.width) && this.width ~= size(im,2))
                im = imresize(im,[this.height,this.width]);
            end
        end
        
        function S = struct(this)
            %STRUCT convert to a struct array
            S = struct(...
                'id',{this.id},...
                'width',{this.width},...
                'height',{this.height},...
                'image_path',{this.image_path},...
                'segmentation',{this.segmentation},...
                'pose',{this.pose});
        end
        
        function obj = copy(this)
            %COPY clone an instance
            obj = cell(size(this));
            for i = 1:numel(this)
                obj{i} = sbu.Photo(this(i).struct);
                if ~isempty(obj{i}.segmentation)
                    obj{i}.segmentation = obj{i}.segmentation.copy;
                end
                if ~isempty(obj{i}.pose)
                    obj{i}.pose = obj{i}.pose.copy;
                end
            end
            obj = reshape([obj{:}],size(this));
        end
        
        function obj = flip(this)
            %FLIP flip left and right of an instance
            if nargout == 0, obj = this; else obj = this.copy; end
            if ~isempty(obj.segmentation)
                obj.segmentation.map = flipdim(obj.segmentation.map,2);
            end
            if ~isempty(obj.pose)
                obj.pose.flip(obj.width);
            end
            obj.image_path = flipdim(obj.image,2);
        end
        
        function obj = resize(this, siz)
            %RESIZE resize an instance
            if nargout == 0, obj = this; else obj = this.copy; end
            if numel(siz)==1 % magnitude specified
                S = [siz,siz];
                siz = round(S .* [obj.height,obj.width]);
            else % size specified
                S = siz ./ [obj.height,obj.width];
            end
            if ~isempty(obj.segmentation)
                obj.segmentation.map = imresize(obj.segmentation.map,siz,'nearest');
            end
            if ~isempty(obj.pose)
                obj.pose.resize(S);
            end
            obj.height = siz(1);
            obj.width = siz(2);
            obj.image_path = imresize(obj.image,siz);
        end
        
        function S = to_uci(this)
            %TO_UCI convert to uci format
            S = struct('im',{this.image_path},...
                'point',cell(1,numel(this)),...
                'labels',cell(1,numel(this)));
            for i = 1:numel(this)
                S(i).point = this(i).pose.to_uci;
                S(i).labels = this(i).segmentation.label_map;
            end
        end
        
        function show(this)
            %SHOW show
            imshow(this.image);
        end
    end
    
    methods (Hidden)
        function obj = saveobj(this)
            %SAVEOBJ callback on save
            obj = this.struct;
        end
    end
    
    methods (Hidden, Static)
        function this = loadobj(obj)
            %LOADOBJ callback on load
            this = sbu.Photo(obj);
        end
    end
end
