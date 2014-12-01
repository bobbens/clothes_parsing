classdef SVM < libsvm.Classifier
    %SVM libsvm classifier
    %
    % Usage:
    %   svm = SVM(labels, samples, 'kernel_type',2,...)
    %   svm = SVM('kernel_type',2,...)
    %   svm = svm.train(labels, samples)
    %   [score, labels, acc] = svm.predict(samples)
    %
    % svm.predict(samples,'noprob') disables probability output even if
    % the model is trained for probability output.
    %
    % Properties:
    % svm_type:    set type of SVM (default 0)
    %     0 -- C-SVC
    %     1 -- nu-SVC
    %     2 -- one-class SVM
    %     3 -- epsilon-SVR
    %     4 -- nu-SVR
    % kernel_type: set type of kernel function (default 2)
    %     0 -- linear: u'*v
    %     1 -- polynomial: (gamma*u'*v + coef0)^degree
    %     2 -- radial basis function: exp(-gamma*|u-v|^2)
    %     3 -- sigmoid: tanh(gamma*u'*v + coef0)
    %     4 -- precomputed kernel (kernel values in instance_matrix)
    % degree:      set degree in kernel function (default 3)
    % gamma:       set gamma in kernel function (default 1/num_features)
    % coef0:       set coef0 in kernel function (default 0)
    % cost:        set the parameter C of C-SVC, epsilon-SVR,
    %              and nu-SVR (default 1)
    % nu:          set the parameter nu of nu-SVC, one-class SVM,
    %              and nu-SVR (default 0.5)
    % epsilon:     set the epsilon in loss function of epsilon-SVR
    %              (default 0.1)
    % cachesize:   set cache memory size in MB (default 100)
    % tolerance:   set tolerance of termination criterion (default 0.001)
    % shrinking:   whether to use the shrinking heuristics, 0 or 1
    %              (default 1)
    % prob:        whether to train a SVC or SVR model for probability
    %              estimates, 0 or 1 (default 0)
    % weight:      set the parameter C of class i to weight*C, for
    %              C-SVC (default 1)
    % n_fold:      n-fold cross validation mode
    % quiet:       quiet mode (no outputs)
    %
    % Packaged by Kota Yamaguchi 2011
    properties
        model       = []
        svm_type    = 0
        kernel_type = 2
        degree
        gamma
        coef0
        cost
        nu
        epsilon
        cachesize
        tolerance
        shrinking
        prob
        weight
        n_fold
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
            if ~isa(y,'double'), y = double(y); end
            this.update_options(varargin{:});
            this.model = libsvm.svmtrain(y, x, this.train_options);
        end
        
        function r = predict(this, y, x, varargin)
            %PREDICT
            if ~isa(y,'double'), y = double(y); end
            r = struct;
            [r.y_hat, r.acc, r.val] = libsvm.svmpredict(y, x,...
                this.model, this.predict_options);
        end
        
        function x = val(this, x)
            %VAL
            r = this.predict(zeros(size(x,1),1),x);
            x = r.val(this.model.Label==1);
        end
        
        function [] = parameter_selection(this, y, x, varargin)
            %PARAMETER_SELECTION search best parameters for C-SVM
            
            % Options
            folds = 5;
            c_grid = 4.^(-1:3);
            g_grid = 2.^(-4:0);
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'folds',  folds = varargin{i+1};
                    case 'cost',   c_grid = varargin{i+1};
                    case 'gamma',  g_grid = varargin{i+1};
                end
            end
            
            if this.svm_type~=0
                warning('LibSVM:parameter_selection',...
                    'parameter selection is not supported for this svm');
            else
                this.n_fold = folds;
                acc = zeros(length(c_grid),length(g_grid));
                for i = 1:length(c_grid)
                    for j = 1:length(g_grid)
                        this.cost  = c_grid(i);
                        this.gamma = g_grid(j);
                        this.train(y,x);
                        acc(i,j) = this.model;
                    end
                end
                [tmp,ind] = max(acc(:));
                [c,g] = ind2sub(size(acc),ind);
                this.cost   = c_grid(c);
                this.gamma  = g_grid(g);
                this.n_fold = [];
                this.model  = [];
            end
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
                        case 'svm_type'
                            c = [c,{sprintf('-s %d',this.svm_type)}];
                        case 'kernel_type'
                            c = [c,{sprintf('-t %d',this.kernel_type)}];
                        case 'degree'
                            c = [c,{sprintf('-d %f',this.degree)}];
                        case 'gamma'
                            c = [c,{sprintf('-g %f',this.gamma)}];
                        case 'coef0'
                            c = [c,{sprintf('-r %f',this.coef0)}];
                        case 'cost'
                            c = [c,{sprintf('-c %f',this.cost)}];
                        case 'nu'
                            c = [c,{sprintf('-n %f',this.nu)}];
                        case 'epsilon'
                            c = [c,{sprintf('-p %f',this.epsilon)}];
                        case 'cachesize'
                            c = [c,{sprintf('-m %d',this.cachesize)}];
                        case 'tolerance'
                            c = [c,{sprintf('-e %f',this.tolerance)}];
                        case 'shrinking'
                            c = [c,{sprintf('-h %d',this.shrinking)}];
                        case 'prob'
                            if this.prob
                                c = [c,{sprintf('-b %d',this.prob)}];
                            end
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
            if ~isempty(this.prob) && this.prob
                c = sprintf('-b %d',this.prob);
            else
                c = '';
            end
        end
    end
    
    methods(Static)
        function this = loadobj(obj)
            %LOADOBJ
            this = libsvm.SVM(obj);
        end
        
        function test_SVM()
            %TEST_SVM
            p = fileparts(mfilename('fullpath'));
            [y,xt] = libsvm.SVM.read([p,filesep,'heart_scale']);
            svm = libsvm.SVM(y,xt);
            r = svm.predict(y, xt);
            disp(r);
        end
    end
    
end

