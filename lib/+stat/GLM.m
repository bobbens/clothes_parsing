classdef GLM < stat.Distribution
    %GLM wrapper class for generalized linear model
    properties(Constant)
        DEFAULT_LINK = struct(...
            'normal',   'identity',...
            'poisson',  'log',...
            'binomial', 'logit',...
            'gamma',    'loglog',...
            'inverse_gaussian', -2 ...
            );
    end
    
    properties
        b
        distr = 'normal'
        link
    end
    
    methods
        function [ this ] = GLM(varargin)
            %GLM constructor
            if numel(varargin)>=2 && isnumeric(varargin{1}) &&...
                    (isnumeric(varargin{2})||islogical(varargin{2}))
                [X,y] = deal(varargin{1},varargin{2});
                varargin = varargin(3:end);
            end
            
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'distr', this.distr = varargin{i+1};
                    case 'link',  this.link = varargin{i+1};
                end
            end
            if isempty(this.link)
                this.link = this.DEFAULT_LINK.(strrep(this.distr,' ','_'));
            end
            ind = find(strcmp(varargin(1:2:numel(varargin)),'distr'));
            if ~isempty(ind), varargin([ind,ind+1]) = []; end
            
            if exist('X','var') && exist('y','var')
                this = this.fit(X,y,varargin{:});
            end
        end
        
        function [ this ] = fit(this, X, y, varargin)
            %FIT
            ind = find(strcmp(varargin(1:2:numel(varargin)),'link'));
            if ~isempty(ind), varargin([ind,ind+1]) = []; end
            
            this.b = glmfit(X,y,this.distr,'link',this.link,varargin{:});
        end
        
        function [ yhat ] = val(this, X, varargin)
            %VAL
            yhat = glmval(this.b,X,this.link,varargin{:});
        end
        
        function [ yhat ] = pdf(this, X, varargin)
            %PDF
            yhat = this.val(X,varargin{:});
        end
    end
    
    methods(Static)
        function [] = test_GLM()
            %TEST_GLM
            import stat.GLM;
            x = [2100 2300 2500 2700 2900 3100 ...
                 3300 3500 3700 3900 4100 4300]';
            y = [ 1  2  0  3  8  8 14 17 19 15 17 21]' ./...
                [48 42 31 34 31 21 23 23 21 16 17 21]';
            glm = GLM(x,y,'distr','binomial');
            yhat = glm.val(x);
            plot(x, y,'o',x, yhat,'-','LineWidth',2);
            fprintf('       x      y   yhat\n');
            fprintf('  %6d %1.4f %1.4f\n', [x,y,yhat]');
            pause;
            close all;
        end
    end
end
