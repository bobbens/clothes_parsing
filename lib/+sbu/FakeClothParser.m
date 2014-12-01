classdef FakeClothParser < sbu.ClothParser
    %FAKECLOTHPARSER fake clothing parser class for baseline comparison
    
    methods
        %% Object initialization
        function [ this ] = FakeClothParser( varargin )
            %FAKECLOTHPARSER create a new clothing parser object
            %
            %    cp = sbu.ClothParser(S)
            %    cp = sbu.ClothParser('PropertyName', propertyValue, ...)
            %
            % ## Input
            % * __S__ struct of property values
            %
            this.initialize(varargin{:});
        end
        
        function [ obj ] = copy(this)
            %COPY deep copy clothing parser object
            %
            %    cloned = cp.copy
            %
            % ## Output
            % * __cloned__ copied clothing parser object
            %
            obj = sbu.FakeClothParser(this.struct);
        end
        
        %% API
        
        function [ Lhat, score, R ] = infer( this, X, varargin )
            %INFER infer labels to the input data sample
            %
            %    [Lhat,score,R] = cp.infer(X, 'ParamName', paramValue, ...)
            %
            % ## Input
            % * __X__ a data sample struct returned by feature_transform
            % * __score__ log probability of the assignment
            % * __R__ output struct returned from libdai
            %
            Lhat = zeros(size(X.L)); % always zero-prediction
            score = -inf;
            R = struct;
        end
        
        function Lhat = infer_each( this, X, varargin )
            %INFER_EACH iterator version of infer method for array input
            %
            % See also sbu.ClothParser.infer
            %
            Lhat = cell(size(X));
            for i = 1:numel(X)
                try
                    Lhat{i} = this.infer(X(i), varargin{:});
                catch e
                    disp(e.getReport);
                end
            end
        end
        
        function train( this, X, varargin )
            %TRAIN learn distributions for input labels and features
            %
            %    cp.train(X, 'ParamName', paramValue, ...)
            %
            % ## Input
            % * __X__ data samples returned by feature_transform
        end
    end
    
    methods (Hidden)
        function [ obj ] = saveobj(this)
            %SAVEOBJ serialize before save
            obj = this.struct;
        end
    end
    
    %%
    methods (Static)
        
        function [ A, CV, Lhat ] = cross_validation( X, varargin )
            %CROSS_VALIDATION apply cross validation to the dataset
            %
            %     [A, R, Lhat] = sbu.ClothParser.cross_validation(X, ...)
            %
            % ## Input
            % * __X__ data samples in a struct array.
            %
            % ## Output
            % * __A__ pixel accuracy of the model
            % * __CV__ cross validation struct containing results and
            %     evaluation
            %
            % ## Parameters
            % * __Nfolds__ number of folds. default 3
            % * __TestSamples__ Optional separate testing samples in the
            %     cross validation. By default, training set is used. Used
            %     in the 'train on true annotation, test on predicted
            %     annotation' scenario
            % * __Verbose__ verbosity flag
            %
            % Also the function takes parameters of the constructor,
            % train, infer
            %
            K = 3;
            verbose = false;
            Xtest = [];
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'Nfolds', K = varargin{i+1};
                    case 'Verbose', verbose = varargin{i+1};
                    case 'TestSamples', Xtest = varargin{i+1};
                end
            end
            
            % Create cross validation structure
            cv_id = mod(1:numel(X),K)+1;
            cv_id = arrayfun(@(x) cv_id==x, 1:K, 'UniformOutput',false);
            CV = struct('ind', cv_id, 'model', cell(1,K),...
                'Lhat', cell(1,K), 'acc', cell(1,K), 'R', cell(1,K));
            
            % Parallel execution
            matlabpool open;
            parfor k = 1:K
                cp = sbu.FakeClothParser(varargin{:});
                cp.train(X(~CV(k).ind), varargin{:});
                CV(k).model = cp.struct;
                if ~isempty(Xtest)
                    x_test = Xtest(CV(k).ind);
                else
                    x_test = X(CV(k).ind);
                end
                Lhat_ = cp.infer_each(x_test, varargin{:});
                L_ = {X(CV(k).ind).L};
                CV(k).Lhat = Lhat_;
                [CV(k).acc, CV(k).R] =...
                    sbu.FakeClothParser.evaluate(L_,Lhat_,{x_test.area});
            end
            matlabpool close;
            
            % Prepare output
            A = mean([CV.acc]);
            Lhat = cell(size(X));
            for k = 1:K, Lhat(CV(k).ind) = CV(k).Lhat; end
            
            if verbose
                for k = 1:K
                    fprintf('Fold %d\n', k);
                    disp(CV(k));
                    disp(CV(k).R);
                end
                fprintf('Average accuracy: %f\n', A);
            end
        end
    end
    
    methods (Static,Hidden)        
        function [ this ] = loadobj(obj)
            %LOADOBJ decode after load
            this = sbu.FakeClothParser(obj);
        end
    end
    
end

