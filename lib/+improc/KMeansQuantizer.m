classdef KMeansQuantizer < improc.Quantizer
    %KMEANSQUANTIZER
    %
    % SYNOPSIS:
    %   q = improc.KMeansQuantizer( ... )
    %   q.learn( vectors )
    %   I = q.quantize( vectors )
    %   H = q.hist( I )
    %   H = q.transform( vectors )
    %
    % KMeansQuantizer assigns scalar group index to a vector with a
    % dictionary trained from training samples.
    % 
    % USAGE:
    %   q = improc.KMeansQuantizer; % construct a new quantizer
    %   q.learn(training_samples);  % learn a dictionary
    %   idx = q.quantize(testing_samples); % assign group indices
    %
    %   h = q.transform(vectors);   % transform vectors into BOW
    %                               % representation
    
    properties
        model     % Learned model
        K = 128   % Number of clusters
        distance = 'euclidean'
        % Distance measure used in the quantizer. It must be one of
        % {'euclidean','cityblock','cosine','correlation','hamming'}
    end
    
    properties(Constant,Hidden)
        DISTANCES = ...
            {'euclidean','cityblock','cosine','correlation','hamming'}
        % Supported distance measures
    end
    
    properties(Dependent,Access=private)
        distance_for_kmeans
    end
    
    methods
        function [ this ] = KMeansQuantizer( varargin )
            %KMEANSQUANTIZER
            this@improc.Quantizer(varargin{:});
        end
        
        function [ obj ] = saveobj(this)
            %SAVEOBJ serialize before save
            obj = this.to_struct;
        end
        
        function [ obj ] = copy(this)
            %COPY
            obj = improc.KMeansQuantizer(this.to_struct);
        end
        
        function [ y ] = transform( this, vectors, varargin )
            %HIST transform set of vectors to a normalized histogram
            idx = this.quantize(vectors, varargin{:});
            y = this.hist(idx);
        end
        
        function [ y ] = hist( this, idx )
            %HIST compute normalized histogram of quantized words
            y = hist(idx(:), 1:this.K);
            n = sum(y);
            if n~=0, y = y/n; end
        end
        
        function [ idx ] = quantize( this, vectors, varargin )
            %QUANTIZE transform vectors to scalar values
            if size(vectors,2)~=size(this.model,2)
                error('KMeansQuantizer:quantize:invalidArg',...
                      'Invalid input vectors');
            end
            
            D = pdist2(this.model,vectors,this.distance);
            [tmp, idx] = min(D,[],1);
            idx = idx(:);
        end
        
        function [ idx ] = learn( this, vectors, varargin )
            %LEARN build a dictionary model from training samples
            
            % Remove invalid option
            ind = find(strcmp(varargin(1:2:length(varargin)),'distance'));
            if ~isempty(ind), varargin([ind,ind+1]) = []; end
            
            % Cluster
            [idx, this.model] = kmeans(vectors, this.K,...
                'distance', this.distance_for_kmeans, varargin{:});
        end
        
        
        % Internal use
        function [] = set.distance(this, d)
            %SET.DISTANCE validates distance name
            if ~any(strcmp(d,this.DISTANCES))
                error('KMeansQuantizer:Distance:Unsupported',...
                      'Unsupported distance measure: %s',d);
            end
            this.distance = d;
        end
        
        function [ name ] = get.distance_for_kmeans(this)
            %GET.DISTANCE_FOR_KMEANS map distance names in pdist to kmeans
            switch this.distance
                case 'euclidean', name = 'sqEuclidean';
                case 'hamming',   name = 'Hamming';
                otherwise,        name = this.distance;
            end
        end
    end
    
    methods(Static)        
        function [] = test_kmeans_quantizer()
            %TEST_KMEANS_QUANTIZER test script
            X = randn(100,3);
            Y = randn(10,3);
            % Train a model
            q = improc.KMeansQuantizer('K',10);
            q.learn(X);
            % Get quantization and histogram
            I = q.quantize(Y);
            I = q.hist(I);
            disp(I);
            % Shorthand
            I = q.transform(Y);
            disp(I);
        end
        
        function [ this ] = loadobj(obj)
            %LOADOBJ decode after load
            this = improc.KMeansQuantizer(obj);
        end
    end
    
end

