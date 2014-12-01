classdef SVM < libsvm.Classifier
    %SVM liblinear classifier
    %
    % Usage:
    %   svm = SVM(labels, samples, 'type',2,...)
    %   svm = SVM('type',2,...)
    %   svm.train(labels, samples)
    %   r = svm.predict(samples)
    %
    % Properties:
    % -s type : set type of solver (default 1)
    % 	0 -- L2-regularized logistic regression (primal)
    % 	1 -- L2-regularized L2-loss support vector classification (dual)
    % 	2 -- L2-regularized L2-loss support vector classification (primal)
    % 	3 -- L2-regularized L1-loss support vector classification (dual)
    % 	4 -- multi-class support vector classification by Crammer and Singer
    % 	5 -- L1-regularized L2-loss support vector classification
    % 	6 -- L1-regularized logistic regression
    % 	7 -- L2-regularized logistic regression (dual)
    % -c cost : set the parameter C (default 1)
    % -e epsilon : set tolerance of termination criterion
    % 	-s 0 and 2
    % 		|f'(w)|_2 <= eps*min(pos,neg)/l*|f'(w0)|_2,
    % 		where f is the primal function and pos/neg are # of
    % 		positive/negative data (default 0.01)
    % 	-s 1, 3, 4 and 7
    % 		Dual maximal violation <= eps; similar to libsvm (default 0.1)
    % 	-s 5 and 6
    % 		|f'(w)|_inf <= eps*min(pos,neg)/l*|f'(w0)|_inf,
    % 		where f is the primal function (default 0.01)
    % -B bias : if bias >= 0, instance x becomes [x; bias]; if < 0, no bias
    %           term added (default -1)
    % -wi weight: weights adjust the parameter C of different classes
    %             (see README for details)
    % -v n: n-fold cross validation mode
    % -q : quiet mode (no outputs)
    %
    % Packaged by Kota Yamaguchi 2011
    properties
        model   = []
        type    = 1
        cost
        epsilon
        bias
        weight
        n_fold
        prob
        quiet
    end
    
    methods
        function this = SVM(varargin)
            %SVM constructor
            this = this@libsvm.Classifier(varargin{:});
        end
        
        function obj = saveobj(this)
            %SAVEOBJ
            obj = this.struct;
        end
        
        function this = train(this, y, x, varargin)
            %TRAIN
            %
            % options:
            %   autoweight: decide
            if ~isa(y,'double'), y = double(y); end
            if ~issparse(x), x = sparse(x); end
            
            % Additional options
            this.update_options(varargin{:});
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'autoweight' % determine weights from data
                        if varargin{i+1}
                            this.weight = this.autoweight(y);
                        end
                end
            end
            
            % Train
            this.model = liblinear.train(y, x, this.train_options);
        end
        
        function r = predict(this, y, x, varargin)
            %PREDICT
            if ~isa(y,'double'), y = double(y); end
            if ~issparse(x), x = sparse(x); end
            r = struct;
            [r.y_hat, r.acc, r.val] = liblinear.predict(y, x,...
                this.model, this.predict_options);
        end
        
        function x = val(this, x)
            %VAL
            this.quiet = true;
            r = this.predict(zeros(size(x,1),1),x);
            x = r.val(this.model.Label==1);
        end
        
        function [] = parameter_selection(this, y, x, varargin)
            %PARAMETER_SELECTION search best parameters for C-SVM
            
            % Options
            folds = 5;
            c_grid = 4.^(-1:3);
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'folds',  folds = varargin{i+1};
                    case 'cost',   c_grid = varargin{i+1};
                end
            end
            
            % CV
            this.n_fold = folds;
            acc = zeros(length(c_grid),1);
            for i = 1:length(c_grid)
                this.cost  = c_grid(i);
                this.train(y,x);
                acc(i) = this.model;
            end
            [tmp,ind] = max(acc(:));
            this.cost   = c_grid(ind);
            this.n_fold = [];
            this.model  = [];
        end
    end
    
    methods (Access = protected)
        function [ c ] = train_options(this)
            %TRAIN_OPTIONS
            c = {};
            props = properties(this);
            for i = 1:numel(props)
                if ~isempty(this.(props{i}));
                    switch props{i}
                        case 'type'
                            c = [c,{sprintf('-s %d',this.type)}];
                        case 'cost'
                            c = [c,{sprintf('-c %f',this.cost)}];
                        case 'epsilon'
                            c = [c,{sprintf('-e %f',this.epsilon)}];
                        case 'bias'
                            c = [c,{sprintf('-B %f',this.bias)}];
                        case 'weight'
                            for i = 1:2:numel(this.weight)
                                c = [c,{sprintf('-w%d %f',...
                                    this.weight(i),this.weight(i+1))}];
                            end
                        case 'n_fold'
                            c = [c,{sprintf('-v %d',this.n_fold)}];
                        case 'quiet'
                            if this.quiet, c = [c,{'-q'}]; end
                    end
                end
            end
            c = [c;repmat({' '},size(c))];
            c = [c{:}];
            if isempty(c), c = ''; end
        end
        
        function [ c ] = predict_options(this)
            %PREDICT_OPTIONS
            c = {};
            if ~isempty(this.prob) && this.prob
                c = [c,{sprintf('-b %d',this.prob)}];
            end
            if this.quiet, c = [c,{'-q'}]; end
            if isempty(c)
                c = '';
            else
                c = [c;repmat({' '},[1,numel(c)])];
                c = [c{1:end-1}];
            end
        end
    end
    
    methods(Static)
        function this = loadobj(obj)
            %LOADOBJ
            this = liblinear.SVM(obj);
        end
        
        function w = autoweight(y)     
            %AUTOWEIGHT return weight vector
            labels = unique(y);
            w = zeros(2,numel(labels));
            for j = 1:numel(labels)
                w(1:2,j) = [labels(j); numel(y)/nnz(y==labels(j))];
            end
            w = w(:)';
        end
        
        function test_SVM()
            %TEST_SVM
            p = fileparts(mfilename('fullpath'));
            [y,xt] = liblinear.SVM.read([p,filesep,'heart_scale']);
            svm = liblinear.SVM(y,xt);
            r = svm.predict(y, xt);
            disp(r);
        end
    end
    
end

