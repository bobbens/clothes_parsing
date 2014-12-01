classdef Distribution
    %DISTRIBUTION mixin
    
    methods
        function [ this ] = Distribution(varargin)
            %DISTRIBUTION default constructor
            if numel(varargin)>=2 && isnumeric(varargin{1}) &&...
                    (isnumeric(varargin{2})||islogical(varargin{2}))
                [X,y] = deal(varargin{1},varargin{2});
                varargin = varargin(3:end);
            end
            if exist('X','var') && exist('y','var')
                this = this.fit(X,y,varargin{:});
            end
        end
        
    end
    
    methods(Abstract)
        this = fit(this, X, y, varargin)
        yhat = val(this, X, varargin)
    end
end
