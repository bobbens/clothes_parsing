classdef Segmentation
    %SEGMENTATION UC Berkeley hierarchical segmentation
    %
    %    [ obj ]    = bsr.Segmentation(I, ...)
    %    [ obj ]    = bsr.Segmentation('ucm',U, ...)
    %    [ obj ]    = bsr.Segmentation(..., 'resize', rsz)
    %
    %    [ bdry ]   = obj.boundary()
    %    [ bdry ]   = obj.boundary(K)
    %
    %    [ labels ] = obj.labels()
    %    [ labels ] = obj.labels(K)
    %
    %    obj.save(filepath)
    %
    %    [ obj ] = bsr.Segmentation.load(filepath)
    %
    %    bsr.make
    %    bsr.make('clean')
    %
    % bsr.Segmentation(im) returns hierarchical image segmentation object
    % obj given an input image I. If you already have computed
    % segmentation, you can supply the hierarchical representation U.
    % After computing segmentation, you can query boundary map with
    % obj.boundary() or segment labels with obj.labels(). These methods
    % take optional threshold value K as an argument. To check what
    % threshold values take an effect, you can use obj.threshold() and get
    % a list of thresholding values from the top. If you choose top N
    % value in this list, you can get N+1 segments in obj.labels().
    %
    % bsr.Segmentation(...) takes optional properties of resizing factor.
    % This parameter takes value between (0,1] and larger values have
    % better segmentation quality but higher computational cost. This
    % resize property can be also set by obj.resize = rsz after the
    % instantiation of an object. The default resize value is 0.5.
    %
    % obj.save() and bsr.Segmentation.load() saves/loads result in/from
    % a mat file respectively. This method is slightly more space
    % efficient than the builtin save() and load() function of matlab.
    %
    % Before using the bsr.Segmentation class, use bsr.compile to build
    % mex functions internally used in bsr package. bsr.compile('clean')
    % will delete these compiled binaries.
    %
    % See also: watershed
    
    % Kota Yamaguchi 2011 <kyamagu@cs.stonybrook.edu>
    
    properties
        ucm             % double-sized boundary map
        resize = 0.5	% resize parameter to reduce computational cost
                        % must be between (0,1]
    end
    
    methods
        function [ this ] = Segmentation(varargin)
            %SEGMENTATION create and initialize an object
            
            % Determine if the initial image is given
            im = [];
            if ~isempty(varargin) && isnumeric(varargin{1})
                im = varargin{1};
                varargin = varargin(2:end);
            end
            
            % Populate properties if given any
            props = properties(this);
            for i = 1:2:length(varargin)
                ind = strcmp(varargin{i},props);
                if any(ind), this.(props{ind}) = varargin{i+1}; end
            end
            
            % Compute globalPb if image is given
            if ~isempty(im)
                this = this.run(im);
            end
        end
        
        function [ this, gpb ] = run(this, im, rsz)
            %RUN compute segmentation of an image 
            if nargin > 2, this.resize = rsz; end
            % compute gpb (note: this takes a minute or so)
            gpb = bsr.globalPb(im, this.resize);
            this.ucm = bsr.contours2ucm(gpb,'doubleSize');
        end
        
        function [ Y ] = prob(this)
            %PROB return boundary probability map
            Y = this.ucm(3:2:end, 3:2:end);
        end
        
        function [ Y ] = boundary(this, k)
            %BOUNDARY return binary boundary map at threshold k
            if nargin < 2, k = 0.05; end
            Y = (this.prob >= k);
        end
        
        function [ Y ] = labels(this, k)
            %LABELS return label matrix at threshold k
            if nargin < 2, k = 0.05; end      
            L = bwlabel(this.ucm < k);
            Y = L(2:2:end, 2:2:end);
        end
        
        function [ Y ] = thresholds(this)
            %THRESHOLDS return possible values of threshold
            Y = this.prob();
            Y = flipud(unique(Y(:)));
        end
        
        function [] = save(this, filepath)
            %SAVE save segmentation object in a mat file
            %  TODO: shall we implement saveobj?
            s = struct('ucm',sparse(this.ucm),'resize',this.resize);
            save(filepath,'-v7.3','-struct','s');
        end
        
        function obj = saveobj(this)
            %SAVEOBJ
            obj = struct('ucm',sparse(this.ucm),'resize',this.resize);
        end
        
        function [] = save_png(this, filepath)
            %SAVE segmentation object in a png file
            imwrite(this.compress_map(this.ucm),filepath,'Gamma',this.resize);
        end
    end
    
    methods(Static)
        function [ this ] = load(filepath)
            %LOAD load segmentation object from a mat file
            %  TODO: shall we implement loadobj?
            s = load(filepath);
            this = bsr.Segmentation('ucm',full(s.ucm),'resize',s.resize);
        end
        
        function [ this ] = loadobj(obj)
            %LOADOBJ
            this = bsr.Segmentation('ucm',full(obj.ucm),'resize',obj.resize);
        end

        function [ this ] = load_png(filepath)
            %LOAD load segmentation object from a png file
            info_ = imfinfo(filepath);
            ucm_ = bsr.Segmentation.uncompress_map(imread(filepath));
            this = bsr.Segmentation('ucm',ucm_,'resize',info_.Gamma);
        end        
        
        function [ Y ] = compress_map( X )
            %COMPRESS_MAP encode
            X = uint32(X*(2^24-1));
            Y = zeros(size(X,1),size(X,2),3,'uint8');
            Y(:,:,1) = bitand(X,uint32(255));
            Y(:,:,2) = bitand(bitshift(X,-8),uint32(255));
            Y(:,:,3) = bitand(bitshift(X,-16),uint32(255));
        end
        
        function [ Y ] = uncompress_map( X )
            %UNCOMPRESS_MAP decode
            X = uint32(X);
            Y = bitor(bitor(X(:,:,1),bitshift(X(:,:,2),8)),bitshift(X(:,:,3),16));
            Y = double(Y)/(2^24-1);
        end
    end
    
end

