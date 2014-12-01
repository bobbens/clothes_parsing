classdef PoseEstimator < handle
   %POSEESTIMATE Pose estimator of [Yang 2011]
   properties
      model
      K = [5 5 5 6 6 6 6 5 5 5 5 5 5 5 5 6 6 6 6 5 5 5 5 5 5 5]
      pa = [0 1 2 3 4 5 6 3 8 9 10 11 12 13 2 15 16 17 18 15 20 21 22 ...
         23 24 25]
      sbin = 4
      name = 'pose_estimator'
      use_labels = false
   end
   
   properties (Constant, Hidden)
      CACHE_PATH = fullfile('tmp','cache')
      MAX_SIZE = [150,150]; % maximum image size [hmax,wmax]
   end
   
   methods
      function this = PoseEstimator(varargin)
         %POSEESTIMATOR construct a new estimator object
         if nargin==1 && ischar(varargin{1})
            varargin{1} = load(varargin{1});
         end
         this.initialize(varargin{:});
      end
      
      function initialize(this, varargin)
         %INITIALIZE populate properties
         if nargin==2 && isstruct(varargin{1})
            S = varargin{1};
         else
            S = struct(varargin{:});
         end
         props = intersect(properties(this),fieldnames(S));
         for i = 1:numel(props)
            this.(props{i}) = S.(props{i});
         end
      end

      function load(this, S)
         this.model  = S.model;
         this.K     = S.K;
         this.pa    = S.pa;
         this.sbin   = S.sbin;
         this.use_labels = S.use_labels;
         this.name   = S.name;
      end
      
      function S = struct(this)
         %STRUCT convert to a struct
         S = struct('model', {this.model},...
                  'K',    {this.K},...
                  'pa',   {this.pa},...
                  'sbin',  {this.sbin},...
                  'use_labels', {this.use_labels},...
                  'name',  {this.name});
      end
      
      function boxes = estimate(this, X, varargin)
         %ESTIMATE estimate a pose for a new image
         %
         %   boxes = estimator.estimate(im, ...)
         %   boxes = estimator.estimate(X, ...)
         %
         % ## Input
         % * __X__ struct array of uci format. See train method
         % * __im__ image
         %
         % ## Output
         % * __boxes__ part bounding boxes
         % 
         % ## Parameters
         % * __labels__ (optional) map of semantic labels associated to
         %    an image in the first argument format
         % * __FixRNG__ flag to fix random number sequence for
         %    reproducibility
         % * __Verbose__ verbosity flag
         %
         labels = [];
         verbose = false;
         fix_rng = true;
         for i = 1:2:numel(varargin)
            switch varargin{i}
               case 'labels', labels = varargin{i+1};
               case 'Verbose', verbose = varargin{i+1};
               case 'FixRNG', fix_rng = varargin{i+1};
            end
         end
         if isnumeric(X)
            X = struct('im', X, 'labels', labels);
         end
         if ~this.use_labels, [X.labels] = deal([]); end
         if fix_rng
            rand('seed',0); % Always start with the same seq
         else
            rand('seed',sum(100*clock));
         end
         
         % resize sample, detect, and resize results
         [X,s] = resize_samples(X, uci.PoseEstimator.MAX_SIZE);
         boxes = testmodel(this.name,this.model,X,'',verbose);
         boxes = resize_boxes(boxes,s);
         if isscalar(boxes), boxes = boxes{:}; end
      end
      
      function train(this, X, varargin)
         %TRAIN train an estimator
         %
         %   estimator.train(X, 'ParamName', paramValue, ...)
         %
         % ## Input
         % * __X__ uci-format struct array that has the following fields
         %    * __im__ path to an image or rgb data
         %    * __points__ 26-by-2 array of body-joint coordinates
         %    * __labels__ (optional) map of semantic labels (clothing)
         %
         % ## Parameters
         % * __FixRNG__ flag to fix random number sequence for
         %    reproducibility
         %
         fix_rng = true;
         for i = 1:2:numel(varargin)
            switch varargin{i}
               case 'UseLabels', this.use_labels = varargin{i+1};
               case 'FixRNG', fix_rng = varargin{i+1};
            end
         end
         
         if ~this.use_labels, [X.labels] = deal([]); end
         if fix_rng
            rand('seed',0); % Always start with the same seq
         else
            rand('seed',sum(100*clock));
         end
         
         % set up data
         cachedir(fullfile(this.CACHE_PATH,this.name));
         X = resize_samples(X, uci.PoseEstimator.MAX_SIZE);
         X = [X,flip(X)];
         train_data = pointtobox(X, this.pa);
         % now train the model
         this.model = trainmodel(this.name,train_data,[],...
            this.K,this.pa,this.sbin);
      end
   end
   
   methods (Static)
      function [PCP, CV, boxes] = cross_validation(X, K, varargin)
         %CROSS_VALIDATION apply cross validation
         %
         %   [pcp, CV, boxes] = uci.PoseEstimator.cross_validation(X, K, ...)
         %
         % ## Input
         % * __X__ struct array in uci sample format
         % * __K__ number of folds in cross validation
         %
         % ## Output
         % * __PCP__ average percentage of correctly located parts
         % * __CV__ cross validation struct containing following info
         %    * __ind__ logical index of testing samples
         %    * __model__ trained model
         %    * __boxes__ detections in box format
         %    * __detRate__ detection rate
         %    * __PCP__ percentage of correctly located parts
         %    * __R__ individual part localization rate. The order is
         %       torso, ul_leg, ur_leg, ll_leg, lr_leg, ul_arm,
         %       ur_arm, ll_arm, lr_arm, head
         % * __boxes__ estimated bounding boxes for each sample
         %
         % ## Options
         % * __UseLabels__ flag to enable context histogram
         % * __CachePrefix__ prefix to the cache path
         % * __TestSamples__ Optional separate testing samples in the
         %    cross validation. By default, training set is used. Used
         %    in the 'train on true annotation, test on predicted
         %    annotation' scenario
         % * __Verbose__ verbosity flag
         %
         
         if nargin < 2, K = 3; end
         
         use_labels_ = false;
         verbose = true;
         Xtest = [];
         cache_prefix = '';
         for i = 1:2:numel(varargin)
            switch varargin{i}
               case 'UseLabels', use_labels_ = varargin{i+1};
               case 'TestSamples', Xtest = varargin{i+1};
               case 'CachePrefix', cache_prefix = varargin{i+1};
               case 'Verbose', verbose = varargin{i+1};
            end
         end
         
         % Create cross validation structure
         cv_id = mod(1:numel(X),K)+1;
         cv_id = arrayfun(@(x) cv_id==x, 1:K, 'UniformOutput',false);
         CV = struct('ind', cv_id, 'model', cell(1,K),'boxes', cell(1,K),...
            'detRate', cell(1,K), ...
            'PCP', cell(1,K), ...
            'R', cell(1,K));
         
         % Parallel execution
         matlabpool open;
         parfor k = 1:K
            name_ = sprintf('%scv%d_k%d_label%d',cache_prefix,K,k,use_labels_);
            pe = uci.PoseEstimator('name', name_, varargin{:});
            pe.train(X(~CV(k).ind), varargin{:});
            CV(k).model = pe.struct;
            pe.model.thresh = min(pe.model.thresh,-2.0);
            if ~isempty(Xtest)
               x_test = Xtest(CV(k).ind);
            else
               x_test = X(CV(k).ind);
            end
            CV(k).boxes = pe.estimate(x_test,varargin{:});
            [CV(k).detRate,CV(k).PCP,CV(k).R] =...
               uci.PoseEstimator.evaluate(X(CV(k).ind),CV(k).boxes);
         end
         matlabpool close;
         PCP = mean([CV.PCP]);
         boxes = cell(size(X));
         for k = 1:K, boxes(CV(k).ind) = CV(k).boxes; end
         
         if verbose
            for k = 1:K, disp(CV(k)); end
            fprintf('Average PCP=%f\n',PCP);
         end
      end

      function [detRate,PCP,R] = evaluate(Xtruth, boxes)
         %EVALUATE compare truth and prediction
         [detRate,PCP,R] = UCI_eval_pcp('', boxes, Xtruth);
      end

      function [pxerr,tags] = evaluate_pxerr( photos, boxes )
         tags = { 'head', 'neck', 'right_shoulder', 'right_elbow', 'right_hand', 'left_shoulder', 'left_elbow', 'left_hand', 'right_hip', 'left_hip', 'right_knee', 'left_knee', 'right_ankle', 'left_ankle' };
         pxerr = zeros(length(photos),length(tags));
         for i=1:length(photos);
            pos    = sbu.Pose( boxes{i}(1,:) );
            pos_gt = photos(i).pose;
            for j=1:length(tags);
               str = tags{j};
               x   = getfield( pos,    [str,'_x'] );
               xgt = getfield( pos_gt, [str,'_x'] );
               y   = getfield( pos,    [str,'_y'] );
               ygt = getfield( pos_gt, [str,'_y'] );
               pxerr(i,j)= norm( [x-xgt,y-ygt] );
            end
         end
      end
      
      function clear_cache
         %CLEAR_CACHE clear training cache
         rmdir([uci.PoseEstimator.CACHE_PATH],'s');
      end
      
      function pe = default_estimator
         %DEFAULT_PARSER load default estimator object trained on
         %Fashionista dataset
         S = load(fullfile(fileparts(mfilename('fullpath')),...
            'default_pose_estimator.mat'));
         pe = S.pose_estimator;
      end
   end
   
   methods (Hidden)
      function obj = saveobj(this)
         %SAVEOBJ callback on save
         obj = this.struct;
      end
   end
   
   methods (Hidden, Static)
      function this = loadobj(obj)
         %LOADOBJ callback on load
         this = uci.PoseEstimator(obj);
      end
   end
end

