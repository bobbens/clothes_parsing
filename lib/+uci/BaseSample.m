classdef BaseSample
    %BASESAMPLE interface of data sample instance for uci.PoseDetector
    properties (Abstract)
        im
        point
        labels
    end
    
    methods
        function initialize(this,varargin)
            %INITIALIZE populate properties
            if nargin==1 && isstruct(varargin{1})
                S = varargin{1};
            else
                S = struct(varargin{:});
            end
            props = intersect(properties(this),fieldnames(S));
            for i = 1:numel(props)
                this.(props{i}) = S.(props{i});
            end
        end
    end
end

