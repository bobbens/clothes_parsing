classdef Pose < handle
    %POSE pose data class
    
    properties (Constant)
        % EDGE definition
        EDGES = struct(...
            'head',          'head',...
            'neck',          'head',...
            'right_shoulder','neck',...
            'right_elbow',   'right_shoulder',...
            'right_hand',    'right_elbow',...
            'left_shoulder', 'neck',...
            'left_elbow',    'left_shoulder',...
            'left_hand',     'left_elbow',...
            'right_hip',     'right_shoulder',...
            'right_knee',    'right_hip',...
            'right_ankle',   'right_knee',...
            'left_hip',      'left_shoulder',...
            'left_knee',     'left_hip',...
            'left_ankle',    'left_knee'...
            )
        % LIMB definition
        LIMBS = struct(...
            'head', struct('s',0.04,'joints',{{'head','neck'}}),...
            'torso', struct('s',0,'joints',{{'neck','right_shoulder','right_hip','left_hip','left_shoulder'}}),...
            'upper_right_arm', struct('s',0.02,'joints',{{'right_shoulder','right_elbow'}}),...
            'upper_left_arm', struct('s',0.02,'joints',{{'left_shoulder','left_elbow'}}),...
            'lower_right_arm', struct('s',0.01,'joints',{{'right_elbow','right_hand'}}),...
            'lower_left_arm', struct('s',0.01,'joints',{{'left_elbow','left_hand'}}),...
            'upper_right_leg', struct('s',0.03,'joints',{{'right_hip','right_knee'}}),...
            'upper_left_leg', struct('s',0.03,'joints',{{'left_hip','left_knee'}}),...
            'lower_right_leg', struct('s',0.02,'joints',{{'right_knee','right_ankle'}}),...
            'lower_left_leg', struct('s',0.02,'joints',{{'left_knee','left_ankle'}})...
            )
    end
    
    properties (Hidden, Constant)
        % Conversion matrix for UCI training data from PARSE annotation
        PARSE_TO_UCI = full(sparse(...
            [1  2  3  4   4   5  6   6   7  8   8   9   9   10 11  11  12 13  13  14 ...
            15 16  16  17 18  18  19 20  20  21  21  22 23  23  24 25  25  26],...
            [14 13 9  9   8   8  8   7   7  9   3   9   3   3  3   2   2  2   1   1 ...
            10 10  11  11 11  12  12 10  4   10  4   4  4   5   5  5   6   6],...
            [1  1  1  1/2 1/2 1  1/2 1/2 1  2/3 1/3 1/3 2/3 1  1/2 1/2 1  1/2 1/2 1 ...
            1  1/2 1/2 1  1/2 1/2 1  2/3 1/3 1/3 2/3 1  1/2 1/2 1  1/2 1/2 1],...
            26,14))
    end
    
    properties % instance variables
        head_x
        head_y
        neck_x
        neck_y
        right_shoulder_x
        right_shoulder_y
        right_elbow_x
        right_elbow_y
        right_hand_x
        right_hand_y
        left_shoulder_x
        left_shoulder_y
        left_elbow_x
        left_elbow_y
        left_hand_x
        left_hand_y
        right_hip_x
        right_hip_y
        left_hip_x
        left_hip_y
        right_knee_x
        right_knee_y
        left_knee_x
        left_knee_y
        right_ankle_x
        right_ankle_y
        left_ankle_x
        left_ankle_y
        score = 0
    end
    
    methods
        function [ this ] = Pose( varargin )
            %SEGMENTATION
            this.initialize(varargin{:});
        end
        
        function initialize( this, varargin )
            %INITIALIZE polymorphic initialization
            %   obj.initialize( filepath )           % load saved pose
            %   obj.initialize( struct )             % import struct
            %   obj.initialize( 'prop', value, ... ) % import list
            
            if isscalar(varargin)
                s = varargin{1};
                
                % Convert import format
                if ischar(s) % Single string. Let's assume a filepath
                    s = load(s);
                elseif isnumeric(s) && isa(s,'uint8') % Image given
                    pe = uci.PoseEstimator.default_estimator;
                    box = pe.estimate(s);
                    s = this.import(box);
                elseif isnumeric(s)
                    s = this.import(s); % Import [Yang11] format
                end
                
                % Copy
                if isstruct(s) && isscalar(s) % Struct: copy fields into properties
                    fields = intersect(fieldnames(s),properties(this));
                    fields = setdiff(fields,{'EDGES','LIMBS'});
                    for i = 1:length(fields)
                        this.(fields{i}) = s.(fields{i});
                    end
                elseif isa(s,'handle') % Object: copy fields into properties
                    fields = intersect(properties(s),properties(this));
                    fields = setdiff(fields,{'EDGES','LIMBS'});
                    for i = 1:length(fields)
                        this.(fields{i}) = s.(fields{i});
                    end
                elseif isstruct(s)
                    for i = 1:numel(s);
                        label = s(i).label;
                        this.([label,'_x']) = s(i).point(1);
                        this.([label,'_y']) = s(i).point(2);
                    end
                end
            else
                % Longer argument list: same with struct
                props = properties(this);
                for i = 1:2:length(varargin)
                    ind = find(strcmp(varargin{i},props),1);
                    if ~isempty(ind)
                        this.(props{ind}) = varargin{i+1};
                    end
                end
            end
        end
        
        function [ x ] = annotation(this)
            %ANNOTATION get compact annotation format
            x = struct(...
                'head',[this.head_x,this.head_y],...
                'neck',[this.neck_x,this.neck_y],...
                'right_shoulder',[this.right_shoulder_x,this.right_shoulder_y],...
                'right_elbow',[this.right_elbow_x,this.right_elbow_y],...
                'right_hand',[this.right_hand_x,this.right_hand_y],...
                'left_shoulder',[this.left_shoulder_x,this.left_shoulder_y],...
                'left_elbow',[this.left_elbow_x,this.left_elbow_y],...
                'left_hand',[this.left_hand_x,this.left_hand_y],...
                'right_hip',[this.right_hip_x,this.right_hip_y],...
                'left_hip',[this.left_hip_x,this.left_hip_y],...
                'right_knee',[this.right_knee_x,this.right_knee_y],...
                'left_knee',[this.left_knee_x,this.left_knee_y],...
                'right_ankle',[this.right_ankle_x,this.right_ankle_y],...
                'left_ankle',[this.left_ankle_x,this.left_ankle_y]...
                );
        end
        
        function [ S ] = struct(this)
            %STRUCT conversion to struct
            props = setdiff(properties(this),{'EDGES','LIMBS'});
            S = struct;
            for i = 1:numel(props), S.(props{i}) = this.(props{i}); end
        end
        
        function obj = copy(this)
            %COPY copy an instance
            obj = sbu.Pose(this.struct);
        end
        
        function [ S ] = to_arr(this)
            %TO_ARR convert to struct array format for feature computation
            s = this.annotation;
            S = struct;
            labels = fieldnames(s);
            for i = 1:length(labels)
                parent = this.EDGES.(labels{i});
                S(i).label = labels{i};
                S(i).point = s.(labels{i});
                S(i).direction = (s.(parent) - S(i).point);
                if all(S(i).direction==0), S(i).direction(2) = -1; end
            end
        end
        
        function I = show(this, varargin)
            %SHOW visualize pose
            I = [];
            col = 'c';
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'image', I = varargin{i+1};
                    case 'color', col = varargin{i+1};
                end
            end
            
            % try find an image if it's a database object
            props = properties(this);
            if isempty(I) && any(strcmp(props,'photo_id'))
                photo = Photo.find(this.photo_id);
                I = photo.image;
            end
            
            % render
            if isempty(I)
                gcf; axis ij;
                xlim([1,400]); ylim([1,600]);
            else
                imshow(I);
            end
            l = this.lines;
            line(l(:,[1,3])',l(:,[2,4])','Color',col,'LineWidth',5);
            if nargout > 0
                drawnow;
                g = getframe(gca);
                I = g.cdata(2:end-1,2:end-1,:);
            end
        end
        
        function l = lines(this)
            %LINES  Output line segments in [x0,y0,x1,y1] row vectors
            childs = fieldnames(this.EDGES)';
            a = this.annotation;
            l = cell(numel(childs),1);
            for i = 1:numel(childs)
                i0 = childs{i};
                i1 = this.EDGES.(childs{i});
                if ~strcmp(i0,i1)
                    l{i} = [a.(i0),a.(i1)];
                end
            end
            l = cat(1,l{:});
        end
        
        function L = limbs(this, s)
            %LIMBS  Output polygons in [x0,...,x3],[y0,...,y3] row vectors
            if nargin < 2
                % Decide width for line segments
                b = this.bbox;
                s = b(4);
            end
            
            a = this.annotation;
            fields = fieldnames(this.LIMBS);
            
            % Create polygons
            X = cell(size(fields));
            Y = cell(size(fields));
            for i = 1:numel(fields)
                name = fields{i};
                S = this.LIMBS.(name);
                XY = cellfun(@(j)a.(j),S.joints,'UniformOutput',false);
                XY = cat(1,XY{:});
                if numel(S.joints)==2
                    % Make a rectangle from a line segment
                    dx = XY(1,1) - XY(2,1);
                    dy = XY(1,2) - XY(2,2);
                    th = cart2pol(dx,dy);
                    mx = - S.s * s .* sin(th); % cos(th+pi/2)
                    my =   S.s * s .* cos(th); % sin(th+pi/2)
                    X{i} = [XY(1,1)+mx,XY(1,1)-mx,XY(2,1)-mx,XY(2,1)+mx]';
                    Y{i} = [XY(1,2)+my,XY(1,2)-my,XY(2,2)-my,XY(2,2)+my]';
                elseif numel(S.joints)>2
                    % Make a polygon
                    X{i} = XY(:,1);
                    Y{i} = XY(:,2);
                end
            end
            L = struct('name',fields,'X',X,'Y',Y);
        end
        
        function [X,Y] = points(this)
            %POINTS Return (X,Y) coordinates of the annotation
            if ~isscalar(this), error('Pose:points','not scalar input'); end
            a = struct2cell(this.annotation);
            X = cellfun(@(c) c(1),a);
            Y = cellfun(@(c) c(2),a);
        end
        
        function P = to_parse(this)
            %TO_PARSE convert to annotation format of PARSE dataset
            P = [...
                this.right_ankle_x,     this.right_ankle_y;...
                this.right_knee_x,      this.right_knee_y;...
                this.right_hip_x,       this.right_hip_y;...
                this.left_hip_x,        this.left_hip_y;...
                this.left_knee_x,       this.left_knee_y;...
                this.left_ankle_x,      this.left_ankle_y;...
                this.right_hand_x,      this.right_hand_y;...
                this.right_elbow_x,     this.right_elbow_y;...
                this.right_shoulder_x,  this.right_shoulder_y;...
                this.left_shoulder_x,   this.left_shoulder_y;...
                this.left_elbow_x,      this.left_elbow_y;...
                this.left_hand_x,       this.left_hand_y;...
                this.neck_x,            this.neck_y;...
                this.head_x,            this.head_y;...
                ];
        end
        
        function P = to_uci(this)
            %TO_UCI convert to uci training format
            P = this.PARSE_TO_UCI * this.to_parse;
        end
        
        function flip(this, width)
            %FLIP flip left and right
            props = properties(this);
            xprops = props(cellfun(@(x)~isempty(x),regexp(props,'.+_x')));
            for i = 1:numel(xprops)
                this.(xprops{i}) = width - this.(xprops{i}) + 1;
            end
            mirrors = {'shoulder','elbow','hand','hip','knee','ankle'};
            xy = {'_x','_y'};
            for i = 1:numel(mirrors)
                for j = 1:numel(xy)
                    lprop = ['left_',mirrors{i},xy{j}];
                    rprop = ['right_',mirrors{i},xy{j}];
                    tmp = this.(lprop);
                    this.(lprop) = this.(rprop);
                    this.(rprop) = tmp;
                end
            end
        end
        
        function resize(this, s)
            %RESIZE resize the annotation
            if isscalar(s), s = [s,s]; end
            props = properties(this);
            xprops = props(cellfun(@(x)~isempty(x),regexp(props,'.+_x')));
            yprops = strrep(xprops,'_x','_y');
            for i = 1:numel(xprops)
                this.(xprops{i}) = s(2) * (this.(xprops{i})-1) + 1;
                this.(yprops{i}) = s(1) * (this.(yprops{i})-1) + 1;
            end
        end
        
        function b = bbox(this, varargin)
            %BBOX Return a bounding box of the model in [x,y,w,h]
            %
            % Options
            %    padding: pixels to pad the box [px,py]
            %        if it's between [0,1], the padding is multiplied by
            %        the width and height of the image size.
            %        if the value is scalar, px = py.
            %    siz: image size. default [w,h] = [400,600]
            %    min_padding: minimum padding size in pixels [px,py]
            %
            
            % Options
            padding = [0.3,0.2];
            min_padding = [55,70];
            siz = [400,600];
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'padding', padding = varargin{i+1};
                    case 'siz', siz = varargin{i+1};
                    case 'min_padding', min_padding = varargin{i+1};
                end
            end
            
            % Get corners
            [X,Y] = this.points;
            x0 = min(X); y0 = min(Y);
            x1 = max(X); y1 = max(Y);
            w = x1 - x0 + 1;
            h = y1 - y0 + 1;
            
            % Adjust padding size
            if isscalar(padding), padding = repmat(padding,1,2); end
            if all(0 < padding & padding < 1)
                padding = round(padding .* [w,h]);
            end
            padding = max(padding,min_padding);
            
            % Get bounding box
            dx = padding(1);
            dy = padding(2);            
            x0 = max(x0-dx,1);
            y0 = max(y0-dy,1);
            x1 = min(x1+dx,siz(1));
            y1 = min(y1+dy,siz(2));
            b = [x0,y0,x1-x0+1,y1-y0+1];
        end
    end
    
    methods (Static)
        function [ this ] = import( box )
            %IMPORT import box data from [yang 2011]
            if isempty(box), error('Pose:import','empty box'); end
            if ischar(box), load(box); end
            d = sbu.Pose.box2skel(box);
            d(1).points = round(d(1).points);
            this = sbu.Pose;
            this.head_x           = d(1).points(1,1);
            this.head_y           = d(1).points(1,2);
            this.neck_x           = d(1).points(2,1);
            this.neck_y           = d(1).points(2,2);
            this.right_shoulder_x = d(1).points(3,1);
            this.right_shoulder_y = d(1).points(3,2);
            this.right_elbow_x    = d(1).points(5,1);
            this.right_elbow_y    = d(1).points(5,2);
            this.right_hand_x     = d(1).points(7,1);
            this.right_hand_y     = d(1).points(7,2);
            this.left_shoulder_x  = d(1).points(15,1);
            this.left_shoulder_y  = d(1).points(15,2);
            this.left_elbow_x     = d(1).points(17,1);
            this.left_elbow_y     = d(1).points(17,2);
            this.left_hand_x      = d(1).points(19,1);
            this.left_hand_y      = d(1).points(19,2);
            this.right_hip_x      = d(1).points(10,1);
            this.right_hip_y      = d(1).points(10,2);
            this.left_hip_x       = d(1).points(22,1);
            this.left_hip_y       = d(1).points(22,2);
            this.right_knee_x     = d(1).points(12,1);
            this.right_knee_y     = d(1).points(12,2);
            this.left_knee_x      = d(1).points(24,1);
            this.left_knee_y      = d(1).points(24,2);
            this.right_ankle_x    = d(1).points(14,1);
            this.right_ankle_y    = d(1).points(14,2);
            this.left_ankle_x     = d(1).points(26,1);
            this.left_ankle_y     = d(1).points(26,2);
        end
        
        function [ d ] = box2skel(boxes)
            %BOX2SKEL Data converter for [yang11] box format
            d = struct('points',cell(1,size(boxes,1)),...
                       'score', cell(1,size(boxes,1)));
            for i = 1:size(boxes,1)
                d(i).points = [mean([boxes(i,1:4:(size(boxes,2)-5));
                                     boxes(i,3:4:(size(boxes,2)-3))]);...
                               mean([boxes(i,2:4:(size(boxes,2)-4));
                                     boxes(i,4:4:(size(boxes,2)-2))])]';
                d(i).score = boxes(i,end);
            end
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
            this = sbu.Pose(obj);
        end
    end
end

