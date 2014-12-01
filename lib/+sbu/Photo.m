classdef Photo < sbu.BasePhoto
    %PHOTO File-based implementation of data instances in clothing parser
    %
    % To construct an sbu.Photo object, give a new image as an argument of
    % the constructor
    %
    %    im = imread('myphoto.jpg');
    %    photo = sbu.Photo(im);
    %
    % By default, sbu.Photo class uses a pose estimator trained on
    % Fashionista dataset to get pose information. To set manually
    % annotated data, create a pose instance first and add it to an
    % argument to the constructor.
    %
    %    pose = sbu.Pose();
    %    pose.head_x = 123;
    %    pose.head_y = 43;
    %    ...
    %    photo = sbu.Photo(im, 'pose', pose);
    %
    % The segmentation object takes the same procedure to set a manual
    % labeling.
    %
    %    S = sbu.Segmentation(im);
    %    save segmentation.mat S
    %
    %    load segmentation.mat S
    %    S.labels = manual_annotation;
    %    photo = sbu.Photo(im, 'segmentation', S);
    %    
    % See also sbu.ClothParser uci.PoseEstimator sbu.Segmentation sbu.Pose
    %
    
    properties
        id = 0
        width     % cached width
        height    % cached height
        image_path = ''
        segmentation
        pose
        metadata = struct % Anything goes here
    end
    
    methods
        function this = Photo(varargin)
            %PHOTO constructs a new object
            if nargin>0;
               this.initialize(varargin{:});
               varargin = varargin(find(cellfun(@ischar,varargin),1):numel(varargin));
               this.process(varargin{:});
            end
        end
        
        function initialize(this, varargin)
            %INITIALIZE populate property values
            
            % Check input args
            if nargin==0, error('sbu:Photo','invalid argument'); end
            I = [];
            if nargin > 0 && isnumeric(varargin{1})
                I = varargin{1};
                varargin = varargin(2:end);
            end
            if isscalar(varargin) && isstruct(varargin{1})
                S = varargin{1};
            else
                c = {varargin{:}};
                S = cell2struct({c{2:2:end}}, {c{1:2:end}}, 2);
                %S = struct(varargin{:});
            end
            if ~isempty(I), S.image_path = I; end
            props = intersect(fieldnames(S),properties(this));
            for i = 1:numel(props)
                this.(props{i}) = S.(props{i});
            end
        end
        
        function process(this, varargin)
            %PROCESS process an image
            if isempty(this.image)
                error('sbu:Photo:invalidArg','image is empty');
            end
            
            verbose = false;
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'Verbose', verbose = varargin{i+1};
                end
            end
            
            if isempty(this.pose)
                if verbose, fprintf('Computing pose...'); tic; end
                this.pose = sbu.Pose(this.image);
                if verbose, fprintf('    %f sec\n',toc); end
            end
            
            if isempty(this.segmentation)
                if verbose, fprintf('Computing superpixels...'); tic; end
                this.segmentation = sbu.Segmentation(this.image);
                if verbose, fprintf('    %f sec\n',toc); end
            end
        end
        
        function x = get.width(this)
            %GET.WIDTH
            if isempty(this.width), this.set_width_and_height; end
            x = this.width;
        end
        
        function x = get.height(this)
            %GET.HEIGHT
            if isempty(this.height), this.set_width_and_height; end
            x = this.height;
        end
    end
    
    methods (Access = private)
        function set_width_and_height(this)
            %SET_WIDTH_AND_HEIGHT
            im = this.image_path;
            if ischar(im) && exist(im,'file'), im = imread(im); end
            this.width = size(im,2);
            this.height = size(im,1);
        end
    end
end
