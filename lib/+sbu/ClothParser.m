classdef ClothParser < handle
   %CLOTHPARSER Clothing parser class
   %
   % The class implements a clothing parser in [Yamaguchi 2012].
   %
   % # Usage
   %
   % ## Prerequisite
   % Build all dependencies with make method:
   %
   %   sbu.make;
   %
   % ## Parsing
   % You need a trained model. The package comes with a default parser
   % trained with Fashionista dataset [Yamaguchi 2012] that you can load
   % with `default_parser` method. If you use a default or already trained
   % a new model, `parse` method will do the job for you:
   %
   %   cp = sbu.ClothParser.default_parser;
   %   clothing_names = {'t-shirt', 'pants'}; % here we want these two to be identified
   %   S = cp.parse(im, clothing_names);
   %
   % `S` is an sbu.Segmentation class object populated with estimated
   % clothing labels. You can access clothing id of superpixels by
   % `labels` property of `S`. The clothing ids can be converted to
   % strings by the `sbu.Fashionista.clothing_names` method:
   %
   %   clothing_ids = S.labels;  % an array of clothing ids associated to
   %                       % each superpixel
   %   clothing_names = sbu.Fashionista.clothing_names(S.labels);
   %
   % You can retrieve mask of each superpixel by `mask_of` method of
   % sbu.Segmentation object:
   %
   %   mask = S.mask_of(1); % return a mask of the first superpixel which
   %                   % is labeled as names{1}
   %
   % Also, you can visualize the result via `show_labels` method in
   % sbu.Segmentation object:
   %
   %   S.show('image', im, 'labels', cp.clothings);
   %
   % To check available clothing items in the parser, you can use
   % `clothing_names` method:
   %
   %   all_clothings = cp.clothing_names;
   % 
   % ## Training
   %
   % You need to prepare an array of sbu.Photo objects with correctly
   % set pose and segmentation objects. The train method takes care of the
   % rest.
   %
   %   photos = sbu.Fashionista.load; % replace this with your data
   %   cp = sbu.ClothParser;
   %   X = sbu.ClothParser.feature_transform(photos);
   %   cp.train(X, 'LabelNames', sbu.Fashionista.clothings);
   %
   % Also cross_validation method measures the performance of the parser.
   % 
   %   [accuracy, CV] = sbu.ClothParser.cross_validation(X);
   %
   % See help sbu.Photo for how to create a new sbu.Photo instance.
   %
   
   properties
      % weights applied to the pairwise potentials
      weights = [0.035020,0.067239]
   end
   
   properties (SetAccess = protected)
      % distribution models for unary potential
      models = struct('label',{},'name',{},'distribution',{})
      % distribution for pairwise potential
      pair_model = struct('prior',[],'distribution',[])
   end
   
   properties (Dependent, SetAccess = protected)
      clothing_names
      clothings
   end
   
   methods
      %% Object initialization
      function [ this ] = ClothParser( varargin )
         %CLOTHPARSER create a new clothing parser object
         %
         %   cp = sbu.ClothParser(S)
         %   cp = sbu.ClothParser('PropertyName', propertyValue, ...)
         %
         % ## Input
         % * __S__ struct of property values
         %
         this.initialize(varargin{:});
      end
      
      function initialize( this, varargin )
         %INITIALIZE initialize clothing parser object with input
         %
         %   cp.initialize(S)
         %   cp.initialize('PropertyName', propertyValue, ...)
         %
         % ## Input
         % * __S__ struct of property values
         % 
         if nargin==2 && isstruct(varargin{1})
            S = varargin{1};
         else
            S = struct(varargin{:});
         end
         props = intersect(properties(this),fieldnames(S));
         for i = 1:numel(props)
            this.(props{i}) = S.(props{i});
         end
         
         % Instantiate class objects
         for j = 1:numel(this.models)
            s = this.models(j).distribution;
            if isstruct(s) && isfield(s,'class')
               this.models(j).distribution = feval(s.class,s);
            end
         end
         s = this.pair_model.distribution;
         if isstruct(s) && isfield(s,'class')
            this.pair_model.distribution = feval(s.class,s);
         end
      end
      
      function [ s ] = struct( this )
         %STRUCT convert sbu.ClothParser to a struct
         %
         %   S = cp.struct
         %   S = struct(cp)
         %
         % ## Output
         % * __S__ struct of property values
         %
         s = struct(...
            'weights',   {this.weights},...
            'models',   {this.models},...
            'pair_model',{this.pair_model}...
            );
         for i = 1:numel(s)
            for j = 1:numel(s(i).models)
               if ~isempty(s(i).models(j).distribution)
                  s(i).models(j).distribution = ...
                     s(i).models(j).distribution.struct;
               end
            end
            if ~isempty(s(i).pair_model.distribution)
               s(i).pair_model.distribution =...
                  s(i).pair_model.distribution.struct;
            end
         end
      end
      
      function [ obj ] = copy(this)
         %COPY deep copy clothing parser object
         %
         %   cloned = cp.copy
         %
         % ## Output
         % * __cloned__ copied clothing parser object
         %
         obj = sbu.ClothParser(this.struct);
      end
      
      function [ names ] = get.clothing_names(this)
         %CLOTHING_NAMES
         names = {this.models.name};
      end
      
      function [ c ] = get.clothings(this)
         %CLOTHINGS
         c = struct('id',   [this.models.label],...
                  'name', {this.models.name});
      end
      
      %% API
      function [ photo, X ] = parse( this, im, clothings, varargin )
         %PARSE parse clothing items in an image
         %
         %   photo = cp.parse(im, clothings, ...)
         %   photo = cp.parse(photo, clothings, ...)
         %   [photo,X] = cp.parse(...)
         %
         % ## Input
         % * __im__ image
         % * __photo__ sbu.Photo object
         % * __clothings__ cell array of clothing names
         %
         % ## Output
         % * __photo__ sbu.Photo object populated with estimated
         %    pose and clothing labels
         % * __X__ processed data sample containing intermediate
         %    results. See feature_transform for detail.
         %
         % ## Parameters
         % * __ID__ Optional id to associate with the output X
         % * __S__ Optional precomputed sbu.Segmentation object
         % * __P__ Optional precomputed sbu.Pose object
         %
         
         verbose = false;
         if nargin < 2, clothings = []; end
         for i = 1:2:numel(varargin)
            switch varargin{i}
               case 'Verbose', verbose = varargin{i+1};
            end
         end
         domain = this.validate_domain(clothings);
         
         % Process pose and segmentation
         photo = sbu.Photo(im, varargin{:});
         S = photo.segmentation;
         
         % Feature computation
         if verbose, fprintf('Computing features...\n'); tic; end
         X = this.feature_transform(photo, varargin{:});
         if verbose, fprintf('   %f sec\n',toc); end
         
         % Inference
         if verbose, fprintf('Infering clothing...\n'); tic; end
         varargin = [{'Domain',domain},varargin];
         S.labels = this.infer(X, varargin{:});
         if verbose, fprintf('   %f sec\n',toc); end
      end
      
      function [ Lhat, score, R ] = infer( this, X, varargin )
         %INFER infer labels to the input data sample
         %
         %   [Lhat,score,R] = cp.infer(X, 'ParamName', paramValue, ...)
         %
         % ## Input
         % * __X__ a data sample struct returned by feature_transform
         % * __score__ log probability of the assignment
         % * __R__ output struct returned from libdai
         %
         % ## Parameters
         % * __Domain__ set of labels to consider in the inference. When
         %    it is `true`, possible labels are determined by the set
         %    of labels contained in the ground truth. When it is
         %    `false`, all known labels are considred (all-way
         %    prediction). If a numeric array is given, only those
         %    labels are considred.
         % * __Verbose__ verbosity flag
         %
         verbose = false;
         domain = true;
         use_trimap = true;
         for i = 1:2:numel(varargin)
            switch varargin{i}
               case 'Verbose', verbose = varargin{i+1};
               case 'Domain', domain = varargin{i+1};
               case 'UseTrimap', use_trimap = varargin{i+1};
            end
         end
         domain = this.decide_domain(domain,unique(X.L)');
         
         if verbose, fprintf('  Constructing factor graph...'); tic; end
         fg = libdai.FactorGraph;
         
         % Set local data
         F = X.F;
         pair_1 = X.P(:,1); pair_2 = X.P(:,2);
         if use_trimap
            trimap = X.trimap;
            F = F(trimap,:);
            ind = trimap(pair_1)&trimap(pair_2);
            ind2 = cumsum(trimap);
            pair_1 = ind2(pair_1(ind));
            pair_2 = ind2(pair_2(ind));
            assert(max(pair_1)<=size(F,1) && max(pair_2)<=size(F,1));
         end
         
         % Set unary factors
         for i = 1:size(F,1)
            p = arrayfun(@(j) this.models(j).distribution.val(F(i,:)),...
               domain)';
            fg.push(i,p);
         end
         
         % Set pairwise factors
         P0 = this.pair_model.prior(domain,domain);
         P0 = P0 ./ sum(P0(:)); % normalize
         P0 = P0.^this.weights(1);
         f = this.pair_feature(F(pair_1,:),F(pair_2,:));
         for i = 1:numel(pair_1)
            p = this.pair_model.distribution.val(f(i,:));
            P = ones(size(P0))*(1-p);
            P(sub2ind(size(P0),1:size(P0,1),1:size(P0,1))) = p;
            fg.push([pair_1(i),pair_2(i)],...
                  P0.*(P.^this.weights(2)));
         end
         if verbose, fprintf('  %f seconds\n',toc); end
         
         % Run
         if verbose, fprintf('  Running inference...'); tic; end
         R = fg.infer('BP');
         if verbose, fprintf('  %f seconds\n',toc); end
         Lhat = arrayfun(@(l) this.models(domain(l)).label, R.qmap);
         score = fg.logP(R.qmap);
         
         % Recover original label vector if trimap is used
         if use_trimap
            Lhat_orig = zeros(numel(X.L),1,class(Lhat));
            Lhat_orig(trimap) = Lhat;
            Lhat = Lhat_orig;
         end

         if verbose
            fprintf('  log probability: %f\n',score);
            fprintf('  segment accuracy: %f\n',nnz(Lhat==X.L)/numel(X.L));
            fprintf('  pixel accuracy: %f\n',sum(X.area(Lhat==X.L))/sum(X.area));
         end
         if all(Lhat==0)
            warning('sbu:ClothParser:infer','all background prediction');
         end
      end
      
      function [Lhat, score] = infer_each( this, X, varargin )
         %INFER_EACH iterator version of infer method for array input
         %
         % See also sbu.ClothParser.infer
         %
         Lhat = cell(size(X));
         score = zeros(size(X));
         for i = 1:numel(X)
            try
               [Lhat{i},score(i)] = this.infer(X(i), varargin{:});
            catch e
               disp(e.getReport);
               Lhat{i} = zeros(size(X(i).L),class(X(i).L));
            end
         end
      end
      
      function train( this, X, varargin )
         %TRAIN learn distributions for input labels and features
         %
         %   cp.train(X, 'ParamName', paramValue, ...)
         %
         % ## Input
         % * __X__ data samples returned by feature_transform
         % 
         % ## Parameters
         % * __Distr__ Distribution function. it is a function hundle
         %    to a constructor of distribution class that takes
         %    (features,labels) pair. By default, it is
         %    `@(x,l) liblinear.SVM(l,x,'prob',true,`
         %    `'type',0,'bias',1,'autoweight',true,'quiet',true)`
         % * __UseTrimap__ flag to enforce null labeling outside of a
         %    bounding box.
         % * __OptimGoal__ Goal of the optimization. default 'none'
         % * __LabelsNames__ Struct array of id and name pairs
         % * __Verbose__ verbosity flag
         %
         verbose = false;
         distr = @(x,l) liblinear.SVM(l,x,'prob',true,...
            'type',0,'bias',1,'autoweight',true,'quiet',true);
         use_trimap = true;
         optim_goal = 'none';
         label_names = [];
         for i = 1:2:numel(varargin)
            switch varargin{i}
               case 'Distr', distr = varargin{i+1};
               case 'Verbose', verbose = varargin{i+1};
               case 'UseTrimap', use_trimap = varargin{i+1};
               case 'OptimGoal', optim_goal = varargin{i+1};
               case 'LabelNames', label_names = varargin{i+1};
            end
         end
         
         % Format data structure
         L = cat(1,X.L);
         F = cat(1,X.F);
         labels = unique(L);
         if isempty(label_names)
            label_names = arrayfun(@(x)num2str(x),labels,...
               'UniformOutput',false);
         elseif isstruct(label_names)
            s = label_names;
            label_names = arrayfun(@(x)s([s.id]==x).name,labels,...
               'UniformOutput',false);
         end
         
         % Drop out-of-bbox segments if using trimap
         if use_trimap
            trimap = cat(1,X.trimap);
            L = L(trimap);
            F = F(trimap,:);
            assert(all(unique(L)==labels));
            for i = 1:numel(X)
               ind = X(i).trimap(X(i).P);
               X(i).P = X(i).P(all(ind,2),:);
            end
            clear trimap ind;
         end
         
         % Learn distributions for unary potential
         if verbose, fprintf('%d segments\n',numel(L)); start = tic; end
         for i = 1:numel(labels)
            ind = [this.models.label]==labels(i);
            if ~any(ind)
               this.models = [this.models,...
                  struct('label',labels(i),...
                        'name', label_names{i},...
                        'distribution',[])];
               ind = [this.models.label]==labels(i);
            end
            if verbose, tic; end
            this.models(ind).distribution = distr(F,L==labels(i));
            if verbose
               fprintf('  label %d: %s (%f seconds)\n',...
                  labels(i),label_names{i},toc);tic;
            end
         end
         
         % Learn pairwise potential
         if verbose, tic; end
         % Compute prior distribution
         m = this.compute_prior(X,labels);
         if verbose, fprintf('prior learning: %f seconds\n',toc); tic; end
         % Train smoothing term
         LP = cell(numel(X),1);
         FP = cell(numel(X),1);
         for i = 1:numel(X)
            LP{i} = X(i).L(X(i).P);
            FP{i} = this.pair_feature(...
               X(i).F(X(i).P(:,1),:),X(i).F(X(i).P(:,2),:));
         end
         LP = cat(1,LP{:});
         FP = cat(1,FP{:});
         s = distr(FP,LP(:,1)==LP(:,2));
         if verbose, fprintf('smoother training: %f seconds\n',toc); end
         this.pair_model = struct('prior',m,'distribution',{s});
         
         if verbose, fprintf('model fitting: %f seconds\n',toc(start)); end
         
         % Learn weights
         if ~strcmp(optim_goal,'none')
            this.weights = fminsearch(@(w)-this.evaluate_model(...
               X, w, optim_goal, varargin{:}),this.weights,...
               optimset('Display','iter'));
         end
      end
      
      function A = evaluate_model( this, X, w, method, varargin )
         %EVALUATE_MODEL evaluate the performance against the data
         this.weights = w;
         [Lhat,score] = this.infer_each(X, varargin{:});
         [acc, R] = sbu.ClothParser.evaluate({X.L},Lhat,{X.area});
         switch method
            case 'acc', A = acc;
            case 'mAR', A = R.mean_avr_rec;
            case 'likelihood', A = sum(score);
            otherwise, A = acc;
         end
      end
      
      function [ P ] = show_potential( this, S, X, varargin )
         %SHOW_POTENTIAL
         %
         %   P = cp.show_potential(S, X)
         %
         % ## Input
         % * __S__ sbu.Segmentation object computed from data sample
         % * __X__ data sample, represented by a struct with the
         %    following fields: 
         %    * __L__ labels
         %    * __F__ features
         %
         % ## Output
         % * __P__ potential maps in M-by-N-by-#(clothings) array where
         %    [M,N] are the size of an image
         %
         
         % Compute potentials
         P = zeros([numel(S.map),numel(this.models)]);
         for i = 1:size(X.F,1)
            m = S.mask_of(i);
            p = arrayfun(@(j) this.models(j).distribution.val(X.F(i,:)),...
               1:numel(this.models));
            P(m(:),:) = repmat(p,[nnz(m(:)),1]);
         end
         P = reshape(P,[size(S.map),1,numel(this.models)]);
         montage(P,'DisplayRange',[]);
      end
      
      function domain = validate_domain(this, domain)
         %VALIDATE_DOMAIN check the validity of the clothing names
         
         % Get all clothings if empty
         all_clothing_names = {this.models.name};
         default_labels = {'null','skin','hair'};
         if nargin < 1 || isempty(domain)
            domain = all_clothing_names;
         end
         % Get id representation of names
         if iscellstr(domain)
            % Ensure default_labels
            ind = cellfun(@(x)~any(strcmp(x,domain)),default_labels);
            domain = [domain,default_labels(ind)];
            
            domain_ = intersect(domain,all_clothing_names);
            if numel(domain_) ~= numel(domain)
               unknown = setdiff(domain,all_clothing_names);
               unknown = sprintf('%s ',unknown{:});
               warning('sbu:ClothParser:validateDomain',...
                  'Unknown clothing names removed: %s',unknown);
            end
            domain = cellfun(@(x)...
               this.models(strcmp(x,{this.models.name})).label,...
               domain_);
         end
         domain = sort(domain(:))';
      end
   end
   
   methods (Hidden)
      function [ obj ] = saveobj(this)
         %SAVEOBJ serialize before save
         obj = this.struct;
      end
   end
   
   %%
   methods (Access=protected)
      function domain = decide_domain(this,domain,labels)
         %DECIDE_DOMAIN parse domain option
         % 
         % true => index_of(labels)
         % false => index_of(any_labels)
         % labels => index_of(labels)
         %
         if islogical(domain)
            if domain
               domain = arrayfun(@(x)find([this.models.label]==x,1),...
                  labels,'UniformOutput',false);
               domain = [domain{cellfun(@(x)~isempty(x),domain)}];
            else
               domain = 1:numel(this.models); % any possible label
            end
         elseif isnumeric(domain)
            domain = arrayfun(@(x)find([this.models.label]==x,1),...
               domain(:)','UniformOutput',false);
            domain = [domain{cellfun(@(x)~isempty(x),domain)}];
         end
      end
   end
   
   %%
   methods (Static)
      function [ X ] = feature_transform( D, varargin )
         %FEATURE_TRANSFORM compute features from data samples
         %
         %   X = sbu.ClothParser.feature_transform(D, ...)
         %
         % ## Input
         % * __D__ sbu.BasePhoto objects
         %
         % ## Output
         % * __X__ struct array of processed data samples. it has the
         %    following fields:
         %    * __L__ labels
         %    * __F__ features
         %    * __P__ segment pair indices; labels => (L(P(i,1)),L(P(i,2)))
         %    * __area__ number of pixels in the segment
         %    * __trimap__ labels indicating {0:bg, 1:probably-fg}
         %
         % ## Parameters
         % * __UsePose__ flag for pose usage. when false, no pose
         %    information is considered in the feature
         % * __Aggregator__ aggregator function used in
         %    improc.ImageProcessor object. By default, it is a
         %    normalized histogram for each computed features.
         % * __Verbose__ verbosity flag
         %
         
         ip = improc.ImageProcessor;
         
         % Options
         nparts = numel(fieldnames(sbu.Pose.EDGES));
         hbins = ip.hist_bins('pose',nparts);
         aggregator = @(x)ip.hist(x,hbins);
         verbose = false;
         use_pose = true;
         for i = 1:2:numel(varargin)
            switch varargin{i}
               case 'Verbose', verbose = varargin{i+1};
               case 'UsePose', use_pose = varargin{i+1};
               case 'Aggregator', aggregator = varargin{i+1};
            end
         end
         
         if verbose
            fprintf('%d images\n',numel(D));
            start = tic; tic;
         end
         
         % Collect features
         X = struct('id',{},'L',{},'F',{},'P',{},'area',{},'trimap',{});
         for i = 1:numel(D)
            im = D(i).image;
            if ischar(im), im = improc.rgbread(im); end
            X(i).id = D(i).id;
            s = D(i).segmentation;
            X(i).L = uint32(s.labels);
            X(i).P = s.pairs;
            X(i).area = s.area;
            
            % Compute feature
            pos = D(i).pose;
            if use_pose
               fmap = ip.process(im, 'pose', pos.to_arr);
            else
               fmap = ip.process(im);
            end
            s.features_from(fmap,aggregator);
            X(i).F = sparse(s.features);
            X(i).trimap = s.is_overlapped(pos.bbox);
            
            if verbose && mod(i, floor(numel(D)/10))==0
               fprintf(' %d images processed (%f seconds)\n',i,toc);
            end
         end
         if verbose
            fprintf('feature computation: %f seconds\n',toc(start));
         end
      end


      
      function [ A, CV, Lhat ] = cross_validation( X, varargin )
         %CROSS_VALIDATION apply cross validation to the dataset
         %
         %    [A, R, Lhat] = sbu.ClothParser.cross_validation(X, ...)
         %
         % ## Input
         % * __X__ data samples in a struct array.
         %
         % ## Output
         % * __A__ pixel accuracy of the model, or if 'UseMAR' is true,
         %    mean average recall of the model.
         % * __CV__ cross validation struct containing results and
         %    evaluation
         %
         % ## Parameters
         % * __Nfolds__ number of folds. default 3
         % * __TestSamples__ Optional separate testing samples in the
         %    cross validation. By default, training set is used. Used
         %    in the 'train on true annotation, test on predicted
         %    annotation' scenario
         % * __Verbose__ verbosity flag
         % * __UseMAR__ flag to use mean average recall for evaluation
         %
         % Also the function takes parameters of the constructor,
         % train, infer
         %
         K = 3;
         verbose = false;
         Xtest = [];
         use_mAR = false;
         for i = 1:2:numel(varargin)
            switch varargin{i}
               case 'Nfolds', K = varargin{i+1};
               case 'Verbose', verbose = varargin{i+1};
               case 'TestSamples', Xtest = varargin{i+1};
               case 'UseMAR', use_mAR = varargin{i+1};
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
            cp = sbu.ClothParser(varargin{:});
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
               sbu.ClothParser.evaluate(L_,Lhat_,{x_test.area});
         end
         matlabpool close;
         
         % Prepare output
         A = mean([CV.acc]);
         if use_mAR
            A = mean(arrayfun(@(r) r.mean_avr_rec, [CV.R]));
         end
         
         Lhat = cell(size(X));
         for k = 1:K, Lhat(CV(k).ind) = CV(k).Lhat; end
         
         if verbose
            for k = 1:K
               fprintf('Fold %d\n', k);
               disp(CV(k));
               disp(CV(k).R);
            end
            fprintf('Average accuracy: %f\n', mean([CV.acc]));
         end
      end
      
      function [A, R] = evaluate(L, Lhat, Weights)
         %EVALUTATE evaluate the structured prediction
         if nargin < 3, Weights = ones(size(L)); end
         
         % remove errorneous predictions
         ind = ~cellfun(@isempty,Lhat);
         L = L(ind);
         Lhat = Lhat(ind);
         Weights = Weights(ind);
         
         % calc stats
         R = stat.clsEval(cat(1,L{:}), cat(1,Lhat{:}),...
            'SampleWeights', cat(1,Weights{:}));
         e = cell(size(L));
         for i = 1:numel(L)
            e{i} = stat.clsEval(L{i},Lhat{i},'SampleWeights',Weights{i});
         end
         R.mean_avr_pre = mean(cellfun(@(x) x.avr_pre,e));
         R.mean_avr_rec = mean(cellfun(@(x) x.avr_rec,e));
         A = R.acc;
      end
      
      function [ R, w, best ] = search_weights(X, varargin)
         %SEARCH_WEIGHTS grid search the best weights of the CRF
         %
         % ## Input
         % * __X__ data samples in a structure array
         %
         % ## Output
         % * __R__ 2D struct array of cross validation results
         % * __w__ weight grids used in the search
         % * __best__ best weights
         %
         % ## Parameters
         % * __wgrid__ grids of weights. default
         %    `{[0,3.^(-9:0)],[0,2.^(-7:1)]}`
         %
         % Also takes parameters of the constructor, cross_validation,
         % train, infer
         %
         % The method internally uses parfor loop.
         %
         
         % option processing
         w = {[0,3.^(-9:0)],[0,2.^(-7:1)]}; % grids
         for i = 1:2:numel(varargin)
            if strcmp(varargin{i},'wgrid'), w = varargin{i+1}; end
         end
         
         % get cross validation accuracy for each grid point
         [w1,w2] = ndgrid(w{1},w{2});
         n = numel(w1);
         A = zeros(size(w1));
         R = cell(size(w1));
         for i = 1:n
            [A(i),R{i}] = sbu.ClothParser.cross_validation(X,...
               'weights',[w1(i),w2(i)],varargin{:});
            fprintf('% 3d / % 3d : w = (%f,%f), A = %f\n',...
               i, n, w1(i), w2(i), A(i));
         end
         R = reshape([R{:}],numel(w{1}),numel(w{2}));
         
         % set the one with the best average pixel accuracy
         [A,ind] = max(A);
         best = [w1(ind),w2(ind)];
         fprintf('Pixel accuracy %f achieved at w = (%f,%f)\n',...
            A, best(1), best(2));
      end
      
      function [wbest,A,flag] = search_weights_by_fminsearch(X, varargin)
         %SEARCH_WEIGHTS_BY_FMINSEARCH search best parameter value of
         %the CRF with fminsearch
         %
         %   [wbest,A,flag] =
         %   sbu.ClothParser.search_weights_by_fminsearch(X)
         %
         % ## Input
         % * __X__ data samples in a structure array
         %
         % ## Output
         % * __wbest__ best weights found
         % * __A__ accuracy at the best weights
         % * __flag__ status of termination
         %
         % ## Parameters
         % * __w0__ initial weights. default [0.035020,0.067239]
         %
         % Also takes parameters of cross_validation, train, infer
         %
         w0 = [0.035020,0.067239];
         for i = 1:2:numel(varargin)
            if strcmp(varargin{i},'w0'), w0 = varargin{i+1}; end
         end
         options = optimset('Display','iter');
         [wbest,A,flag] = fminsearch(@(x)-sbu.ClothParser.cross_validation(X,...
            'weights',x,varargin{:}),w0,options);
      end
      
      function f = pair_feature(f1, f2)
         %PAIR_FEATURE compute a pairwise feature
         %
         %   f = sbu.ClothParser.pair_feature(f1,f2)
         %
         % ## Input
         % * __f1__, __f2__ unary features
         %
         % ## Output
         % * __f__ pairwise feature
         %
         f = [0.5*(f1+f2),abs(f1-f2)];
      end
      
      function H = compute_prior(X, labels)
         %COMPUTE_PRIOR obtain the prior matrix
         %
         %   H = sbu.ClothParser.compute_prior(X, labels)
         %
         % ## Input
         % * __X__ data samples in a struct array
         % * __labels__ vector of known labels in a numeric array
         %
         % ## Output
         % * __H__ #labels-by-#labels array of prior potential
         %
         N = zeros(numel(labels));
         H = zeros(numel(labels));
         for i = 1:numel(X)
            L = X(i).L(X(i).P); % label pairs
            I = arrayfun(@(l)find(l==labels,1),L); % convert to index
            ind = I(:,1)>I(:,2);
            I(ind,:) = fliplr(I(ind,:)); % ascending order
            h = accumarray(I,1,[numel(labels),numel(labels)]);
            h = h ./ numel(ind); % get the ratio
            assert(all(h(:)<=1));
            H = H + h;
            N = N + double(h>0);
         end
         H = H ./ N; % get the average ratio
         H(isnan(H)) = min(H(~isnan(H))) / 2; % avoid zero prob
         H(H<H') = H((H<H')'); % make it symmetric
      end
      
      function cp = default_parser
         %DEFAULT_PARSER load default parser object trained on
         %Fashionista dataset
         S = load(fullfile(fileparts(mfilename('fullpath')),...
            'default_clothing_parser.mat'));
         cp = S.cp;
      end
   end
   
   methods (Static,Hidden)
      function [ this ] = loadobj(obj)
         %LOADOBJ decode after load
         this = sbu.ClothParser(obj);
      end
   end
   
end

