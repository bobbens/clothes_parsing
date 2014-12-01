classdef GMM < stat.Distribution
    %GMM gaussian mixture model
    properties
        obj
        K   = 4
        reg = 0.2   % regularization coefficients when singular
    end
    
    methods
        function [ this ] = GMM(varargin)
            %GMM constructor
            if nargin>0 && isnumeric(varargin{1})
                X = varargin{1};
                varargin = varargin(2:end);
                if islogical(varargin{1})
                    X = X(varargin{1},:);
                    varargin = varargin(2:end);
                end
            end
            
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'K', this.K = varargin{i+1};
                    case 'reg', this.reg = varargin{i+1};
                end
            end
            
            if exist('X','var')
                this = this.fit(X,varargin{:});
            end
        end
        
        function [ this ] = fit(this, X, varargin)
            %FIT
            if size(X,1)<size(X,2)
                warning('GMM:TooFewN',...
                    'X must have more rows than columns.');
                this.obj = stat.MVN(X, varargin{:});
            else
                I = kmeans(X,this.K);
                this.obj = gmdistribution.fit(X,this.K,...
                    'Regularize',this.reg,varargin{:},...
                    'Start',I);
            end
        end
        
        function [ p ] = val(this, X, varargin)
            %VAL
            p = this.obj.pdf(X);
        end
        
        function [ p ] = pdf(this, X, varargin)
            %PDF
            p = this.val(X,varargin{:});
        end
    end
    
    methods(Static)
        function [] = test_GMM()
            %TEST_GMM
            import stat.GMM;
        end
    end
end
