classdef MVN < stat.Distribution
    %MVN multivariate normal distribution
    properties
        mu
        inv_sigma
        det_sigma
        reg = 1   % regularization coefficients when singular
    end
    
    properties(Dependent)
        sigma
    end
    
    methods
        function [ this ] = MVN(varargin)
            %GLM constructor
            if nargin>0 && isnumeric(varargin{1})
                X = varargin{1};
                varargin = varargin(2:end);
                if numel(varargin)>0 && islogical(varargin{1})
                    X = X(varargin{1},:);
                    varargin = varargin(2:end);
                end
            end
            
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'mu', this.mu = varargin{i+1};
                    case 'inv_sigma',  this.inv_sigma = varargin{i+1};
                    case 'det_sigma',  this.det_sigma = varargin{i+1};
                    case 'sigma', this.sigma = varargin{i+1};
                end
            end
            
            if exist('X','var')
                this = this.fit(X,varargin{:});
            end
        end
        
        function this = set.sigma(this, S)
            %SET.SIGMA
            if isvector(S), S = diag(S); end
            this.inv_sigma = inv(S);
            this.det_sigma = det(S);
        end
        
        function [ this ] = fit(this, X, varargin)
            %FIT
            if nargin > 2 && strcmp(varargin{1},'diag')
                this.sigma = var(X,0,1);
            else
                s = cov(X);
                if rcond(s)<0.1, s = s + this.reg*eye(size(X,2)); end
                this.sigma = s;
            end
            this.mu = mean(X);
        end
        
        function [ p ] = val(this, X, varargin)
            %VAL
            X = X - repmat(this.mu,[size(X,1),1]);
            p = exp(-0.5*sum(X*this.inv_sigma.*X,2)) /...
                sqrt((2*pi)^(size(X,2))*this.det_sigma);
        end
        
        function [ p ] = pdf(this, X, varargin)
            %PDF
            p = this.val(X,varargin{:});
        end
    end
    
    methods(Static)
        function [] = test_MVN()
            %TEST_MVN
            import stat.MVN;
            mvn = MVN(randn(100,2));
            [X,Y] = meshgrid(-3:.2:3,-3:.2:3);
            P = mvn.pdf([X(:),Y(:)]);
            P = reshape(P,size(X));
            surf(P);
        end
    end
end
