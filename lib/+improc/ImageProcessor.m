classdef ImageProcessor
    %IMAGEPROCESSOR take an image and compute feature map
    %
    % USAGE:
    %   im = imread('/path/to/image.png');
    %   ip = ImageProcessor;
    %   fmap = ip.process(im);
    %
    %   fmap: M-by-N-by-D array of pixelwise dense feature
    %
    % Kota Yamaguchi
    
    properties
        gf = improc.GaborFilter
    end
    
    methods
        function [ this ] = ImageProcessor( varargin )
            %IMAGEPROCESSOR
            props = properties(this);
            for i = 1:2:length(varargin)
                ind = strcmp(varargin{i},props);
                if any(ind), this.(props{ind}) = varargin{i+1}; end
            end
        end
        
        function [ Y ] = process( this, X, varargin )
            %PROCESS preprocess an image for feature representation
            import improc.ColorSpace;
            if ischar(X), X = imread(X); end
            Y = struct;
            Y.rgb = im2double(X);
            Y.lab = ColorSpace.rgb2lab(X) / 108.883;
            Y.texture = this.gf.apply(X);
            [Y.x,Y.y] = meshgrid((1:size(X,2))/size(X,2),...
                                 (1:size(X,1))/size(X,1));
            
            for i = 1:2:length(varargin)
                if strcmp(varargin{i},'pose')
                    p = varargin{i+1};
                    x = zeros(size(X,1),size(X,2),numel(p));
                    y = zeros(size(X,1),size(X,2),numel(p));
                    for j = 1:numel(p)
                        [x(:,:,j),y(:,:,j)] = ...
                            this.xymap([size(X,2),size(X,1)],...
                                       p(j).point, p(j).direction);
                    end
                    Y.pose = cat(3,x,y);
                end
            end
            
            % Pack result in an array
            Y = struct2cell(Y);
            Y = cat(3,Y{:});
        end
        
        function [ e ] = hist_bins(this, varargin)
            %HIST_BINS
            ngrid = 11;
            pose = 14;
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'pose', pose = varargin{i+1};
                    case 'ngrid', ngrid = varargin{i+1};
                end
            end
            
            this.gf.size;
            e = [...
                repmat({linspace(0,1,ngrid)},1,3),...               % rgb
                repmat({linspace(0,1,ngrid)},1,3),...               % lab
                repmat({linspace(-.5,.5,ngrid)},1,this.gf.size),... % gf
                repmat({linspace(-1,1,ngrid)},1,2),...              % xy
                repmat({linspace(-1,1,ngrid)},1,2*pose)...          % pose
                ];
        end
        
        function [ x ] = hist(this, x, bins)
            %HIST
            if nargin < 3, bins = this.hist_bins; end
            c = cell(1,size(x,2));
            for i = 1:size(x,2)
                c{i} = hist(x(:,i),bins{i});
            end
            x = [c{:}];
            n = sum(x(:));
            if n~=0, x = x / n; end % normalize
        end
    end
    
    methods(Static)
        function [ x, y ] = xymap( d, t, direction )
            %XYMAP xy coordinates relative to point t, orientation theta
            %          d: [width, height]
            %          t: [x0, y0]
            %  direction: [dx, dy]
            
            % rotate 90 degree so that vertical axis corresponds to the
            % direction of the parent
            theta = - pi/2 + atan2(direction(2)/d(2),direction(1)/d(1));
            % here, the origin is the lower left corner of an image
            [x,y] = meshgrid((1:d(1))/d(1),(1:d(2))/d(2));
            Z = [x(:),y(:),ones(numel(x),1)] * ...  % original coordinates
                [1,0,0;0,1,0;-t(1)/d(1),-t(2)/d(2),1] *...   % translation
                [ cos(theta), sin(theta),0;...      % rotation
                 -sin(theta), cos(theta),0;...
                           0,          0,1];
            x(:) = Z(:,1);
            y(:) = Z(:,2);
        end
    end
    
end

