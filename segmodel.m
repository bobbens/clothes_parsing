classdef segmodel
   %EXPERIMENT Main experimental protocol
   
   properties (Constant, Hidden)
      TMP = 'tmp' % path to the temp directory
      FEATURES_FILE = 'features.mat';
      FEATURES_POSE_FILE = 'features_pose.mat';
      SHAPELETS_FILE = 'shapelets_features.mat';
      SHAPELETS_POSE_FILE = 'shapelets_features_pose.mat';
      ZFEATURES_FILE = 'Zfeatures.mat';
      CPMC_FILE = 'cpmc.mat';
      CPMC_UNARIES_54_FILE = 'cpmc_unaries_54.mat';
      CPMC_UNARIES_29_FILE = 'cpmc_unaries_29.mat';
      CPMC_UNARIES_11_FILE = 'cpmc_unaries_11.mat';
      TAMARA_UNARIES_54_FILE = 'tamara_unaries_54.mat';
      TAMARA_UNARIES_29_FILE = 'tamara_unaries_29.mat';
      TAMARA_UNARIES_11_FILE = 'tamara_unaries_11.mat';
      TAMARA_UNARIES_POSE_54_FILE = 'tamara_unaries_pose_54.mat';
      TAMARA_UNARIES_POSE_29_FILE = 'tamara_unaries_pose_29.mat';
      TAMARA_UNARIES_POSE_11_FILE = 'tamara_unaries_pose_11.mat';
      LIMBS_FILE = 'limbs_e5_w40.mat';
      JOINTS_FILE = 'joints_s2.5.mat'
      LABELS_FILE = 'labels_limbs.mat';
      POSE_FILE = 'pose.mat';
      SIMILARITY_FEATURES_FILE = 'photos_similarity.mat';
      SIMILARITY_FILE = 'seg_similarity.mat';
      SIMILARITY_POSE_FILE = 'seg_similarity_pose.mat';
   end

   properties
      % Paths
      PROFILE = '0.16';
      SEGMENTPATH = 'poseseg/output/oracle_segments/SuperSegment0.16/';
      SEGMENTS_FILE = '%06d_oracle.mat';
      FEATPATH  = 'poseseg/features/0.16/';
      FEATURES = [];
      FEATURES_POSE = [];
      SEGMENTS = [];
      SHAPELETS = [];
      SHAPELETS_POSE = [];
      ZFEATURES = [];
      CPMCFEATS = '/share/data/poseseg/output/cpmccls29/MyMeasurements/';
      CPMC     = [];
      CPMC_UNARIES_54 = [];
      CPMC_UNARIES_29 = [];
      CPMC_UNARIES_11 = [];
      TAMARA_UNARIES_54 = [];
      TAMARA_UNARIES_29 = [];
      TAMARA_UNARIES_11 = [];
      TAMARA_UNARIES_POSE_54 = [];
      TAMARA_UNARIES_POSE_29 = [];
      TAMARA_UNARIES_POSE_11 = [];
      LIMBS = [];
      JOINTS = [];
      LABELS = [];
      POSE = [];
      SIMILARITY_FEATURES = [];
      SIMILARITY = [];
      SIMILARITY_POSE = [];
      SIMILARITIES = 'data/wn_%s.mat';
      GTMAPS   = 'poseseg/output/segments_groundtruth/%06d.mat';
      CPMC_PRED = 'poseseg/output/cpmccls/preds/%06d.mat';
      CPMC_APPROX = 'poseseg/output/cpmccls/MySegmentsMat/CPMC_segms_150_sp_approx/%06d.mat';
      % Distribution parameters
      NUM_MACHINES   = 5; % Number of machines to use

      % Data
      use_fashionista_v2 = false; % What dataset to use

      % Feature parameters
      pa    = [0 1 2 3 4 5 6 3 8 9 10 11 12 13 2 15 16 17 18 15 20 21 22 23 24 25];
      shape_size = repmat( 1.2, 1, 56 );

      ngrid = 11; % Number of histogram bins for features
      symm_limbs = [  3  4 15 16;
                      4  5 16 17;
                      5  6 17 18;
                      6  7 18 19;
                      3  8 15 20;
                      8  9 20 21;
                      9 10 21 22;
                     10 11 22 23;
                     11 12 23 24;
                     12 13 24 25;
                     13 14 25 26 ];
      symm_joints = [ 3:14;15:26 ]';
      symm_limbs_tree = [ 1  2;
                          2  3;
                          3  4;
                          1  5;
                          5  6;
                          6  7;
                          7  8;
                          8  9;
                          9 10;
                         10 11 ];

      % Indexes and features
      photos_truth   = []; % Ground truth annotated photos
      clothings      = []; % Clothing labels
      clothings11    = [];
      clothings29    = [];
      clothings54    = [];
      map_54_target  = [];
      map_29_target  = [];
      map_11_target  = [];
      trainids       = []; % Indexes of photos to use to train
      testids        = []; % Indexes of photos to use to test
      X              = []; % Segment unary features

      % Internal state data
      class_occurences = []; % Average times a class appears per image
      weights        = []; % Weights for MRF
      wlabels        = {}; % Labels for each weight
      Params         = []; % Parameters for dSP

      % Dataset information
      class_group    = 'fine'; % How to group, can be 'fine', 'coarse', 'medium' or 'none'/[]

      % MRF Parameters
      class_disimilarity = []; % Class disimilarity
      use_disimilarity = 'lch'; % Disimilarity to use: wup, lch or path
      use_disim_loss = true; % Use disimilarity loss
      use_class_loss = true; % Weight loss by class
      use_pixel_loss = false; % Whether to use pixel area loss or just 1/0 for segments.
      loss_background = 16; % In times the largest disimilarity metric
      loss_superpixels = 10;
      % Generally doesn't improve results at cost of more weights
      use_fullbias   = false; % Whether or not to use bias terms for all classes or just background
      % Very important
      use_logreg_54  = false; % Whether to use logistic regression or feature vector directly
      use_logreg_54_max = false;
      use_logreg_29  = true;
      use_logreg_11  = false;
      use_logreg_11_neg = true;
      % Very important
      use_shapelets  = true; % Whether or not to use shapelets
      % Improves about 0.5% over baseline
      use_cpmc       = true; % Whether or not to use CMPC as unary
      use_cpmc_01    = false; % Whether or not to use CMPC as just a 0-1 value (not weighted by score),
                              % Seems to be about the same as not weighing
      use_cpmc_logreg_54 = false; % Use the logistic regression for cpmc unaries
      use_cpmc_logreg_54_max = false;
      use_cpmc_logreg_29 = true;
      use_cpmc_logreg_11 = false;
      use_cpmc_logreg_11_neg = true;
      % Similarity pairwise connections
      use_similarity = true; % Use similarity
      similarity_threshold = 0.85; % Percent threshold
      similarity_full = true; % Per class weights
      similarity_mintree = true; % Use minimum spanning tree
      % Symmetry connections 
      use_limbs      = true; % Use limbs
      use_limbs_symm = true; % Use symmetry potentials across different parts
      use_limbs_symm_limbs = true; % Whether to use limbs or joints
      use_limbs_pose = true;
      use_limbs_pose_pos = true;
      use_limbs_pose_neg = true;
      use_limbs_merge = true; % Merge and don't keep separate states
      symm_extension = 5;%3; % Gaussian length multiplier
      symm_width     = 40; % Gaussian width (pixels)
      symm_thresh    = 0.01;%0.2; % Threshold range for superpixels to be considered
      symm_sigma     = 2.5; % Gaussian width (scaling) for the joint
      loss_L         = 1e1; % Loss on L

      % Real pose
      % We have to recalculate both shapelets and tamara unaries for this
      use_real_pose = false; % Evaluate on test using real pose

      % Pot multipliers
      pot_seg_logreg54 = 1e3;
      pot_seg_logreg29 = 1e3;
      pot_seg_logreg11 = 1e3;
      pot_seg_cpmc_logreg54 = [1e3, 1e3, 1e3];
      pot_seg_cpmc_logreg29 = [1e3, 1e3, 1e3];
      pot_seg_cpmc_logreg11 = [1e3, 1e3, 1e3];
      pot_seg_bias   = 1e1;
      pot_seg_shapelets = 1e3;
      pot_seg_cpmc   = 1e3;
      pot_similarity = 1e2; % Similarity pot
      pot_L_bias     = 1e1;
      pot_L_L        = 1e1;
      pot_L_L_pose   = 1e1;
      pot_L_L_pose_neg = 1e1;
      pot_seg_L      = 1e1;

      % Superpixels
      use_gt_superpixels = false; % Whether or not to use tho annotated ground truth superpixels
      %seg_threshold  = 0.16; % Threshold for superpixels
   end
   
   methods (Static)
      function RMRF = main()
         e = segmodel

         % Z unaries
         e = e.train_misc_unaries();

         % MRF
         e = e.train_MRF();
         RMRF = e.test_MRF_segmentation();
       end

      % Exports the MRF so it can be run with mdSP
      function [e, Pred] = MRF_import( path )
         if exist( [path,'/model.mat'], 'file' );
            data  = load( [path,'/model.mat'], 'model' );
            e     = segmodel( data.model );
         else
            e     = segmodel;
            warning('No model file found');
         end
         if exist( [path,'/weights'], 'file' );
            w        = ReadOutput( [path,'/weights'] );
            e.weights = w;
         end
         if exist( [path,'/predict'], 'file' );
            [w,Pred] = ReadOutput( [path,'/predict'] );
         else
            Pred     = [];
         end
      end

      % Low memory version
      function [R] = MRF_evaluate( path )
         Rpath = [path,'/results_seg.mat'];
         if ~exist(Rpath,'file');
            [e, Pred] = segmodel.MRF_import( path );
         end
         if exist(Rpath,'file');
            load( Rpath, 'R' );
         else
            R = e.test_MRF_segmentation( Pred );
            save( Rpath, 'R' );
         end

         fid = fopen( [path,'/results'], 'w' );
         fprintf( fid, 'acc:     %.4f\n',   R.acc );
         fprintf( fid, 'avr_pre: %.4f\n',   R.avr_pre );
         fprintf( fid, 'avr_rec: %.4f\n',   R.avr_rec );
         fprintf( fid, 'avr_voc: %.4f\n\n', R.avr_voc );
         fclose(  fid );
      end
   end

   % Most of the MRF internal building routines are here
   % The features are in the 0-1 range unless reweighted
   % The loss for class occurences is roughly in the range 1-1000
   % The loss for class disimilarity is in the range 1-2.5 (wup)
   methods (Access = protected)
      % Builds the unary nodes for the segments within the MRF
      function [S, numvars, IdsSeg, wlabels] = build_MRF_segments( this, S, O, seg, numvars, X, wlabels )
         L  = zeros(size(seg.labels));
         for i=1:numel(this.clothings);
            L( ismember(seg.labels, this.clothings(i).id) ) = i;
         end
         assert(~any(L==0))
         A  = X.area; % Area of each segment in pixels

         % Unaries for each segment
         IdsSeg = 1:O.nsegments; % Vector to map segment to region id
         for j=1:O.nsegments;
            S.Regions{j}.c_r              = 1; % why?
            % Each unary only uses the segment variable
            S.Regions{j}.VariableIndices  = [j];
            % Higher order factors and connections would go here
            S.Regions{j}.Parents          = [];
            % Features
            r  = 0; % Number of variables (assumed unique)
            fc = []; % Column values
            fr = []; % Row values
            fv = []; % Values
            % Using Logistic Regression as main feature
            if this.use_logreg_29;
               map  = this.map_29_target;
               F  = (map .* repmat(X.tamara_unaries_29(j,:)',1,size(map,2)))';
               iscol = (sum(F,1)>0);
               ncol = sum(iscol);
               [a,b,c] = find( F(:,iscol)  );
               fc = [fc a'];
               fr = [fr r+b'];
               fv = [fv this.pot_seg_logreg29.*c'];
               r  = r+ncol;
               if j==1;
                  for nc=1:size(F,2);
                     if ~iscol(nc); continue; end
                     wlabels = {wlabels{:} sprintf('log reg 29 class=%s', this.clothings29(nc).name)};
                  end
               end
            end
            if this.use_logreg_11;
               map  = this.map_11_target;
               F  = (map .* repmat(X.tamara_unaries_11(j,:)',1,size(map,2)))';
               iscol = (sum(F,1)>0);
               ncol = sum(iscol);
               [a,b,c] = find( F(:,iscol)  );
               fc = [fc a'];
               fr = [fr r+b'];
               if this.use_logreg_11_neg;
                  fv = [fv -this.pot_seg_logreg11.*(1-c)'];
               else
                  fv = [fv this.pot_seg_logreg11.*c'];
               end
               r  = r+ncol;
               if j==1;
                  for nc=1:size(F,2);
                     if ~iscol(nc); continue; end
                     wlabels = {wlabels{:} sprintf('log reg 11 class=%s', this.clothings11(nc).name)};
                  end
               end
            end
            if this.use_logreg_54;
               map  = this.map_54_target;
               F  = (map .* repmat(X.tamara_unaries_54(j,:)',1,size(map,2)))';
               nF = size(F,2);
               iscol = (sum(F,1)>0);
               F  = F(:,iscol);
               if this.use_logreg_54_max;
                  dF = diag(max(F,[],2));
                  dF = dF( :, sum(dF>0,1)>0 );
                  [a,b,c] = find( dF );
               else
                  [a,b,c] = find( F );
               end
               fc = [fc a'];
               fr = [fr r+b'];
               fv = [fv this.pot_seg_logreg54.*c'];
               r  = r+numel(c);
               if j==1;
                  for nc=1:nF;
                     if ~iscol(nc); continue; end
                     wlabels = {wlabels{:} sprintf('log reg 54 class=%s', this.clothings54(nc).name)};
                  end
               end
            end
            % Using cpmc logistic regressions
            if this.use_cpmc_logreg_29;
               for k=1:3;
                  if this.pot_seg_cpmc_logreg29(k)==0;
                     continue;
                  end
                  map  = this.map_29_target;
                  F  = (map .* repmat(X.cpmc_unaries_29{k}(j,:)',1,size(map,2)))';
                  iscol = (sum(F,1)>0);
                  ncol = sum(iscol);
                  [a,b,c] = find( F(:,iscol)  );
                  fc = [fc a'];
                  fr = [fr r+b'];
                  fv = [fv this.pot_seg_cpmc_logreg29(k).*c'];
                  r  = r+ncol;
                  if j==1;
                     for nc=1:size(F,2);
                        if ~iscol(nc); continue; end
                        wlabels = {wlabels{:} sprintf('CPMC log reg 29 feat=%d, class=%s', k, this.clothings29(nc).name )};
                     end
                  end
               end
            end
            if this.use_cpmc_logreg_11;
               for k=1:3;
                  if this.pot_seg_cpmc_logreg11(k)==0;
                     continue;
                  end
                  map  = this.map_11_target;
                  F  = (map .* repmat(X.cpmc_unaries_11{k}(j,:)',1,size(map,2)))';
                  iscol = (sum(F,1)>0);
                  ncol = sum(iscol);
                  [a,b,c] = find( F(:,iscol)  );
                  fc = [fc a'];
                  fr = [fr r+b'];
                  if this.use_cpmc_logreg_11_neg
                     fv = [fv -this.pot_seg_cpmc_logreg11(k).*(1-c)'];
                  else
                     fv = [fv this.pot_seg_cpmc_logreg11(k).*c'];
                  end
                  r  = r+ncol;
                  if j==1;
                     for nc=1:size(F,2);
                        if ~iscol(nc); continue; end
                        wlabels = {wlabels{:} sprintf('CPMC log reg 11 feat=%d, class=%s', k, this.clothings11(nc).name )};
                     end
                  end
               end
            end
            if this.use_cpmc_logreg_54;
               for k=1:3;
                  if this.pot_seg_cpmc_logreg54(k)==0;
                     continue;
                  end
                  map  = this.map_54_target;
                  F  = (map .* repmat(X.cpmc_unaries_54{k}(j,:)',1,size(map,2)))';
                  nF = size(F,2);
                  iscol = (sum(F,1)>0);
                  F  = F(:,iscol);
                  if this.use_logreg_54_max;
                     dF = diag(max(F,[],2));
                     dF = dF( :, sum(dF>0,1)>0 );
                     [a,b,c] = find( dF );
                  else
                     [a,b,c] = find( F );
                  end
                  fc = [fc a'];
                  fr = [fr r+b'];
                  fv = [fv this.pot_seg_cpmc_logreg54(k).*c'];
                  r  = r+numel(c);
                  if j==1;
                     for nc=1:nF;
                        if ~iscol(nc); continue; end
                        wlabels = {wlabels{:} sprintf('CPMC log reg 54 feat=%d, class=%s', k, this.clothings54(nc).name )};
                     end
                  end
               end
            end
            % Using bias terms for ALL the classes
            if this.use_fullbias;
               fc = [fc 1:O.NUM_CLASSES];
               fr = [fr r+(1:O.NUM_CLASSES)];
               fv = [fv this.pot_seg_bias.*ones(1,O.NUM_CLASSES)];
               r  = r+O.NUM_CLASSES;
               if j==1;
                  for nc=1:O.NUM_CLASSES;
                     wlabels = {wlabels{:} sprintf('bias class=%s', this.clothings(nc).name )};
                  end
               end
            % Only bias term for the background
            else
               fc = [fc 1];
               fr = [fr r+1];
               fv = [fv this.pot_seg_bias];
               r  = r+1;
               if j==1;
                  wlabels = {wlabels{:} 'background bias'};
               end
            end
            % Using shapelets as an additional feature
            if this.use_shapelets;
               % Add to sparse matrix
               fc = [fc 1:O.NUM_CLASSES];
               fr = [fr r+(1:O.NUM_CLASSES)];
               fv = [fv this.pot_seg_shapelets.*X.shapelets_score(j,:)];
               r  = r+O.NUM_CLASSES;
               if j==1;
                  for nc=1:O.NUM_CLASSES;
                     wlabels = {wlabels{:} sprintf('shapelet class=%s', this.clothings(nc).name )};
                  end
               end
            end
            if this.use_cpmc;
               cpmc = X.cpmc_score(j,:);
               if this.use_cpmc_01;
                  cpmc = cpmc ./ sum(cpmc);
               end
               fc = [fc 1:O.NUM_CLASSES];
               fr = [fr repmat(r+1, 1, O.NUM_CLASSES)];
               fv = [fv this.pot_seg_cpmc.*[cpmc(2) repmat( cpmc(1), 1, O.NUM_CLASSES-1)]];
               r  = r+1;
               if j==1;
                  wlabels = {wlabels{:} 'cpmc'};
               end
            end
            % Build the features as a sparse matrix
            S.Regions{j}.r    = numvars + 1:r;
            %fc,fr,fv,O.NUM_CLASSES,r
            %[numel(fc),numel(fr),numel(fv),O.NUM_CLASSES,r]
            S.Regions{j}.pot  = sparse( fc, fr, fv, O.NUM_CLASSES, r );

            % Loss function - VERY IMPORTANT
            if this.use_disim_loss;
               S.Regions{j}.Loss          = this.class_disimilarity(L(j), :); % based on class
            else
               S.Regions{j}.Loss          = ones( O.NUM_CLASSES, 1 ); % 0/1 Loss
            end
            if this.use_pixel_loss;
               S.Regions{j}.Loss          = S.Regions{j}.Loss .* A(j) ./ sum(A); % Loss is pixel area normalized
            end
            % Penalize more uncommon classes gotten wrong
            if this.use_class_loss;
               max_occ  = max(this.class_occurences);
               seg_occ  = this.class_occurences(L(j));
               loss_fac = max_occ ./ seg_occ;
               S.Regions{j}.Loss    = S.Regions{j}.Loss .* loss_fac;
               %S.Regions{j}.Loss    = S.Regions{j}.Loss .* (max(max_occ/1e6, this.class_occurences) ./ max_occ)';
            end
            S.Regions{j}.Loss       = S.Regions{j}.Loss .* this.loss_superpixels;
            S.Regions{j}.Loss(L(j)) = 0; % True annotation
         end
         numvars = numvars + numel(S.Regions{1}.r);
         % Observations
         S.Observation(1:O.nsegments)    = L;
      end

      function [S, numvars, wlabels] = build_MRF_similarity( this, S, O, seg, numvars, X, IdsSeg, wlabels )
         ss = X.similarity_score( (X.similarity_score(:,3) > this.similarity_threshold), : );
         if this.similarity_mintree;
            smin      = ss;
            smin(:,3) = 1-smin(:,3);
            if numel(ss)>0;
               nMST  = grMinSpanTree( smin );
               ss    = ss(nMST,:);
            end
         end
         NC = O.NUM_CLASSES;
         diag_sing = sparse( find(diag(ones(NC,1))>0), ones(NC,1), ones(NC,1), NC*NC, 1 );
         diag_full = sparse( 1:(NC+1):(NC*NC), 1:NC, ones(NC,1), NC*NC, NC );
         n = numel(S.Regions);
         for i=1:size(ss,1);
            n = n+1;
            idSeg1 = IdsSeg(ss(i,1));
            idSeg2 = IdsSeg(ss(i,2));
            S.Regions{n}.c_r     = 1;
            S.Regions{n}.VariableIndices = [idSeg1 idSeg2];
            if this.similarity_full;
               S.Regions{n}.r       = numvars + (1:NC);
               S.Regions{n}.pot     = (this.pot_similarity*ss(i,3)) .* diag_full;
            else
               S.Regions{n}.r       = numvars + 1;
               S.Regions{n}.pot     = (this.pot_similarity*ss(i,3)) .* diag_sing;
            end
            S.Regions{n}.Loss    = sparse( NC*NC, 1 );
            S.Regions{n}.Parents = [];
            S.Regions{idSeg1}.Parents = [S.Regions{idSeg1}.Parents n];
            S.Regions{idSeg2}.Parents = [S.Regions{idSeg2}.Parents n];
            if i==1;
               if this.similarity_full;
                  for j=1:O.NUM_CLASSES;
                     wlabels = {wlabels{:} sprintf( 'Similarity class=%s', this.clothings(j).name)};
                  end
               else
                  wlabels = {wlabels{:} 'Similarity'};
               end
            end
         end
         if this.similarity_full;
            numvars = numvars + NC;
         else
            numvars = numvars + 1;
         end
      end

      function [S, numvars, wlabels] = build_MRF_seg_symmetry_merge( this, S, O, seg, numvars, X, IdsSeg, IdsL, wlabels, L_cooc )

         score = X.L_score;
         limbs = this.symm_limbs;
         limbIds = zeros( size(limbs,1), 1 );
         %wlabels = {wlabels{:} 'Limb bias'};
         diag_pot = sparse( find(diag(ones(O.NUM_CLASSES,1))>0), ...
               ones(O.NUM_CLASSES,1), ones(O.NUM_CLASSES,1), ...
               O.NUM_CLASSES*O.NUM_CLASSES, 1 );

         %numvars = numvars+1;
         for i=1:O.NUM_CLASSES;
            wlabels = {wlabels{:} sprintf( 'Limb bias=%d', i ) };
         end
         n = numel(S.Regions);
         for i=1:size(limbs,1);
            % Create the limb
            n = n+1;
            limbIds(i,1) = n;
            limbIds(i,2) = n-1;
            S.Regions{n}.c_r       = 1;
            S.Regions{n}.VariableIndices   = [ IdsL(i,1) ];
            S.Regions{n}.r         = numvars + (1:O.NUM_CLASSES);
            S.Regions{n}.pot       = this.pot_L_bias .* diag(ones( O.NUM_CLASSES, 1 ));
            %S.Regions{n}.pot       = zeros( O.NUM_CLASSES, 1 );
            %S.Regions{n}.pot(1)    = this.pot_L_bias;
            S.Regions{n}.Parents   = [];

            obs = this.X(i).L_labels(i,3);
            S.Observation( IdsL(i,1) ) = obs;
            S.Regions{n}.Loss        = this.loss_L .* ones( 1, O.NUM_CLASSES );
            S.Regions{n}.Loss( obs ) = 0;
         end
         numvars = numvars+O.NUM_CLASSES;

         % Here superpixels are connected to the auxiliary meta-limb variables
         n = numel(S.Regions);
         nbase = n;
         for i = 1:size(limbs,1);
            inrange1 = find( score(:,i,1) > this.symm_thresh )';
            inrange2 = find( score(:,i,2) > this.symm_thresh )';
            wlabels = {wlabels{:} sprintf('Limb Pair=%d',i)};
            for i1=inrange1;
               idSeg1 = IdsSeg(i1);
               n     = n+1;
               S.Regions{n}.c_r  = 1;
               S.Regions{n}.VariableIndices = [IdsL(i,1) idSeg1];
               S.Regions{n}.r    = numvars+i;
               S.Regions{n}.pot  = score(i1,i,1).*this.pot_seg_L.*diag_pot;
               S.Regions{n}.pot(1,1) = 0;
               % No lossVV
               S.Regions{n}.Loss = sparse( O.NUM_CLASSES*O.NUM_CLASSES, 1 );
               S.Regions{n}.Parents = [];
               % Connect as parent to other regions
               S.Regions{limbIds(i,1)}.Parents = [S.Regions{limbIds(i,1)}.Parents n];
               S.Regions{idSeg1}.Parents       = [S.Regions{idSeg1}.Parents n];
            end
            for i2=inrange2;
               idSeg2 = IdsSeg(i2);
               n     = n+1;
               S.Regions{n}.c_r  = 1;
               S.Regions{n}.VariableIndices = [IdsL(i,1) idSeg2];
               S.Regions{n}.r    = numvars+i;
               S.Regions{n}.pot  = score(i2,i,2).*this.pot_seg_L.*diag_pot;
               S.Regions{n}.pot(1,1) = 0;
               % No lossVV
               S.Regions{n}.Loss = sparse( O.NUM_CLASSES*O.NUM_CLASSES, 1 );
               S.Regions{n}.Parents = [];
               % Connect as parent to other regions
               S.Regions{limbIds(i,1)}.Parents = [S.Regions{limbIds(i,1)}.Parents n];
               S.Regions{idSeg2}.Parents       = [S.Regions{idSeg2}.Parents n];
            end
         end
         numvars = numvars + size(limbs,1);

         % Here we set up pose structure
         if this.use_limbs_pose;
            n = numel(S.Regions);
            for i = 1:size(this.symm_limbs_tree)
               n  = n+1;

               id1 = IdsL( this.symm_limbs_tree(i,1), 1 );
               id2 = IdsL( this.symm_limbs_tree(i,2), 1 );
               i1 = limbIds( this.symm_limbs_tree(i,1), 1 );
               i2 = limbIds( this.symm_limbs_tree(i,2), 1 );
               S.Regions{n}.c_r  = 1;
               S.Regions{n}.VariableIndices = [id1 id2];
               C = sum(L_cooc(:,:,i,:),4);
               S.Regions{n}.r    = [];
               S.Regions{n}.pot  = [];
               if this.use_limbs_pose_neg;
                  wlabels = {wlabels{:} sprintf('Limb Pose Negative=%d',i)};
                  numvars = numvars+1;
                  S.Regions{n}.r    = [S.Regions{n}.r numvars];
                  S.Regions{n}.pot  = [S.Regions{n}.pot, this.pot_L_L_pose_neg .* ...
                        (-reshape( C==0, O.NUM_CLASSES*O.NUM_CLASSES, 1 ))];
               end
               if this.use_limbs_pose_pos;
                  wlabels = {wlabels{:} sprintf('Limb Pose Positive=%d',i)};
                  R = 1./sum(C);
                  R(isinf(R))=0;
                  numvars = numvars+1;
                  S.Regions{n}.r    = [S.Regions{n}.r numvars];
                  S.Regions{n}.pot  = [S.Regions{n}.pot, this.pot_L_L_pose .* ...
                        reshape( C .* repmat( R, size(C,1), 1 ), O.NUM_CLASSES*O.NUM_CLASSES, 1 )];
               end
               S.Regions{n}.Loss = sparse( O.NUM_CLASSES*O.NUM_CLASSES, 1 );
               S.Regions{n}.Parents = [];
               S.Regions{i1}.Parents = [S.Regions{i1}.Parents n];
               S.Regions{i2}.Parents = [S.Regions{i2}.Parents n];
            end
         end
      end
   end

   methods
      % Helper
      function this = segmodel(varargin)
         %EXPERIMENT Constructor for experiment object
         if nargin==1 && isstruct(varargin{1})
            S = varargin{1};
         else
            S = struct(varargin{:});
         end
         props = intersect(properties(this),fieldnames(S));
         for i = 1:numel(props)
            this.(props{i}) = S.(props{i});
         end

         if strcmpi( this.PROFILE, '0.16' );
            this.SEGMENTPATH = 'poseseg/output/oracle_segments/SuperSegment0.16/';
            this.FEATPATH    = 'poseseg/features/0.16/';
            this.SEGMENTS   = [this.SEGMENTPATH this.SEGMENTS_FILE];
            this.CPMCFEATS = '/share/data/poseseg/output/cpmccls29_0.16/MyMeasurements/';
         elseif strcmpi( this.PROFILE, '0.10' );
            this.SEGMENTPATH = 'poseseg/output/oracle_segments/SuperSegment0.10/';
            this.FEATPATH    = 'poseseg/features/0.10/';
            this.SEGMENTS    = [this.SEGMENTPATH this.SEGMENTS_FILE];
            this.CPMCFEATS = '/share/data/poseseg/output/cpmccls29/MyMeasurements/';
         elseif strcmpi( this.PROFILE, 'gt' );
            this.use_gt_superpixels = true;
            this.FEATPATH    = 'poseseg/features/gt/';
            this.FEATURES    = 'tmp/ttsplit29/clothing_parsing_features.mat';
            this.SEGMENTS    = [];
            this.CPMCFEATS = '/share/data/poseseg/output/cpmccls29gt/MyMeasurements/';
         elseif strcmpi( this.PROFILE, '0.16-11' );
            this.SEGMENTPATH = 'poseseg/output/oracle_segments/SuperSegment0.16/';
            this.FEATPATH    = 'poseseg/features/0.16_11c/';
            this.SEGMENTS    = [this.SEGMENTPATH this.SEGMENTS_FILE];
            this.CPMCFEATS = '/share/data/poseseg/output/cpmccls29_0.16/MyMeasurements/';
         elseif strcmpi( this.PROFILE, '0.10-11' );
            this.SEGMENTPATH = 'poseseg/output/oracle_segments/SuperSegment0.10/';
            this.FEATPATH    = 'poseseg/features/0.10_11c/';
            this.SEGMENTS    = [this.SEGMENTPATH this.SEGMENTS_FILE];
            this.CPMCFEATS = '/share/data/poseseg/output/cpmccls29/MyMeasurements/';
         elseif strcmpi( this.PROFILE, 'gt-11' );
            this.use_logreg_29      = false;
            this.use_logreg_11      = true;
            this.use_cpmc_logreg_29 = false;
            this.use_cpmc_logreg_11 = true;
            this.use_gt_superpixels = true;
            this.class_group = 'medium';
            this.FEATPATH  = 'poseseg/features/gt_11c/';
            this.SEGMENTS  = [];
            this.CPMCFEATS = '/share/data/poseseg/output/cpmccls29gt/MyMeasurements/';
         elseif strcmpi( this.PROFILE, 'iccv' )
            this.use_fashionista_v2 = true;
            this.use_logreg_54      = true;
            this.use_logreg_29      = false;
            this.use_logreg_11      = false;
            this.use_cpmc_logreg_54 = true;
            this.use_cpmc_logreg_29 = false;
            this.use_cpmc_logreg_11 = false;
            this.use_gt_superpixels = true;
            this.class_group = 'none';
            this.use_disimilarity = 'lch_v2';
            this.FEATPATH  = 'poseseg/features/iccv/';
            this.SEGMENTS  = [];
            this.CPMC_PRED = 'poseseg/output/cpmccls_v2.0_fg/preds/%06d.mat';
            this.CPMC_APPROX = 'poseseg/output/cpmccls_v2.0_fg/MySegmentsMat/CPMC_segms_150_sp_approx/%06d.mat';
            this.CPMCFEATS = '/share/data/poseseg/output/cpmccls54_v2.0/SP_gt/MyMeasurements/';
            this.GTMAPS   = 'poseseg/output/segments_groundtruth_v2/%06d.mat';
         elseif strcmpi( this.PROFILE, 'iccv-0.16' )
            this.use_fashionista_v2 = true;
            this.use_logreg_54      = true;
            this.use_logreg_29      = false;
            this.use_logreg_11      = false;
            this.use_cpmc_logreg_54 = true;
            this.use_cpmc_logreg_29 = false;
            this.use_cpmc_logreg_11 = false;
            this.class_group = 'none';
            this.use_disimilarity = 'lch_v2';
            this.SEGMENTPATH = 'poseseg/output/oracle_segments/SuperSegment0.16_v2/';
            this.FEATPATH    = 'poseseg/features/iccv-0.16/';
            this.SEGMENTS    = [this.SEGMENTPATH this.SEGMENTS_FILE];
            this.CPMC_PRED = 'poseseg/output/cpmccls_v2.0_fg/preds/%06d.mat';
            this.CPMC_APPROX = 'poseseg/output/cpmccls_v2.0_fg/MySegmentsMat/CPMC_segms_150_sp_approx/%06d.mat';
            this.CPMCFEATS = '/share/data/poseseg/output/cpmccls54_v2.0/SP_0.16/MyMeasurements/';
            this.GTMAPS   = 'poseseg/output/segments_groundtruth_v2/%06d.mat';
         elseif strcmpi( this.PROFILE, 'iccv-0.10' )
            this.use_fashionista_v2 = true;
            this.use_logreg_54      = true;
            this.use_logreg_29      = false;
            this.use_logreg_11      = false;
            this.use_cpmc_logreg_54 = true;
            this.use_cpmc_logreg_29 = false;
            this.use_cpmc_logreg_11 = false;
            this.class_group = 'none';
            this.use_disimilarity = 'lch_v2';
            this.SEGMENTPATH = 'poseseg/output/oracle_segments/SuperSegment0.10_v2/';
            this.FEATPATH    = 'poseseg/features/iccv-0.10/';
            this.SEGMENTS    = [this.SEGMENTPATH this.SEGMENTS_FILE];
            this.CPMC_PRED = 'poseseg/output/cpmccls_v2.0_fg/preds/%06d.mat';
            this.CPMC_APPROX = 'poseseg/output/cpmccls_v2.0_fg/MySegmentsMat/CPMC_segms_150_sp_approx/%06d.mat';
            this.CPMCFEATS = '/share/data/poseseg/output/cpmccls54_v2.0/SP_0.10/MyMeasurements/';
            this.GTMAPS   = 'poseseg/output/segments_groundtruth_v2/%06d.mat';
         else
            error( ['Unknown profile ' this.PROFILE] )
         end

         % Do again to override stuff that could be profile-based set
         props = intersect(properties(this),fieldnames(S));
         for i = 1:numel(props)
            this.(props{i}) = S.(props{i});
         end

         % Set paths
         this.FEATURES = [this.FEATPATH this.FEATURES_FILE];
         this.FEATURES_POSE = [this.FEATPATH this.FEATURES_POSE_FILE];
         this.SHAPELETS = [this.FEATPATH this.SHAPELETS_FILE];
         this.SHAPELETS_POSE = [this.FEATPATH this.SHAPELETS_POSE_FILE];
         this.ZFEATURES = [this.FEATPATH this.ZFEATURES_FILE];
         this.CPMC = [this.FEATPATH this.CPMC_FILE];
         this.CPMC_UNARIES_54 = [this.FEATPATH this.CPMC_UNARIES_54_FILE];
         this.CPMC_UNARIES_29 = [this.FEATPATH this.CPMC_UNARIES_29_FILE];
         this.CPMC_UNARIES_11 = [this.FEATPATH this.CPMC_UNARIES_11_FILE];
         this.TAMARA_UNARIES_54 = [this.FEATPATH this.TAMARA_UNARIES_54_FILE];
         this.TAMARA_UNARIES_29 = [this.FEATPATH this.TAMARA_UNARIES_29_FILE];
         this.TAMARA_UNARIES_11 = [this.FEATPATH this.TAMARA_UNARIES_11_FILE];
         this.TAMARA_UNARIES_POSE_54 = [this.FEATPATH this.TAMARA_UNARIES_POSE_54_FILE];
         this.TAMARA_UNARIES_POSE_29 = [this.FEATPATH this.TAMARA_UNARIES_POSE_29_FILE];
         this.TAMARA_UNARIES_POSE_11 = [this.FEATPATH this.TAMARA_UNARIES_POSE_11_FILE];
         this.LIMBS = [this.FEATPATH this.LIMBS_FILE];
         this.JOINTS = [this.FEATPATH this.JOINTS_FILE];
         this.LABELS = [this.FEATPATH this.LABELS_FILE];
         this.POSE = [this.FEATPATH this.POSE_FILE];
         this.SIMILARITY_FEATURES = [this.FEATPATH this.SIMILARITY_FEATURES_FILE];
         this.SIMILARITY = [this.FEATPATH this.SIMILARITY_FILE];
         this.SIMILARITY_POSE = [this.FEATPATH this.SIMILARITY_POSE_FILE];

         % Load features and stuff
         this = this.features_load();

         % how many Armijo iterations to perform.
         Params.ArmijoIterations = 50;
         % regularization parameter C used in the equations
         Params.C                = 1;
         % whether to erase the messages after every iteration. Useful if you
         % want to learn CRFs or structured Support Vector Machines the
         % traditional way.
         Params.CRFEraseMessages = 0; % Default 1, Alex says 0
         % how many CRF iterations during one CCCP iteration
         Params.CRFIterations    = 500; % Defaults 50
         % how many message passing iterations during one CRF iteration
         Params.CRFMPIterations  = 1; % Default 50, Alex says 1
         Params.MPIterations     = 500; % Default 50
         % how many message passing iterations during one CRF iteration
         Params.epsilon          = 1;
         % regularization parameter p used in the equations
         Params.p                = 2;
         % whether to reuse the previously computed messages during the Armijo
         % iterations or whether to go with the traditional approach, i.e.,
         % perform inference again
         Params.ReuseMessagesForF = 1; % Default 0, Alex says 1
         % FULL VERBOSITY
         Params.Verbosity        = 0;

         this.Params = Params;
      end

      function [s] = struct( this )
         props = properties(this);
         s = struct;
         for i=1:numel(props);
            s.(props{i}) = this.(props{i});
         end
      end

      function [this] = features_load( this )
         % Load photos
         if this.use_fashionista_v2;
            this.photos_truth = sbu.Fashionista.load_v2;

            % Calculate split
            load( 'data/split_v2.mat' );
            this.trainids  = trainids;
            this.testids   = testids;
         else
            this.photos_truth = sbu.Fashionista.load;
         
            % Calculate split
            load( 'data/split.mat' );
            this.trainids  = trainids;
            this.testids   = testids;
         end
         if isstr( this.SEGMENTS );
            for i=1:numel(this.photos_truth);
               load( sprintf(this.SEGMENTS, i), 'labels', 'map' );
               segstr = struct( 'labels', labels, 'map', map );
               this.photos_truth(i).segmentation = sbu.Segmentation( segstr );
            end
         end

         % Use real pose
         if this.use_real_pose;
            data = load( this.POSE );
            for i = 1:numel( this.photos_truth )
               this.photos_truth(i).pose = sbu.Pose(data.boxes{i});
            end
         end

         % Load features
         if ~exist( this.FEATURES, 'file' );
            warning( ['Features file "',this.FEATURES,'" is missing!'] );
         else
            F     = load( this.FEATURES, 'Xtruth' );
            this.X = rmfield( F.Xtruth, 'F' ); % Get rid of featurers
            % Tamara's logistic regression
            if this.use_real_pose;
               TUNARIES_54 = this.TAMARA_UNARIES_POSE_54;
               TUNARIES_29 = this.TAMARA_UNARIES_POSE_29;
               TUNARIES_11 = this.TAMARA_UNARIES_POSE_11;
            else
               TUNARIES_54 = this.TAMARA_UNARIES_54;
               TUNARIES_29 = this.TAMARA_UNARIES_29;
               TUNARIES_11 = this.TAMARA_UNARIES_11;
            end
            if this.use_logreg_54;
               data = load( TUNARIES_54, 'score' );
               for i=1:numel(this.X);
                  this.X(i).tamara_unaries_54 = data.score{i};
               end
            end
            if this.use_logreg_29;
               data = load( TUNARIES_29, 'score' );
               for i=1:numel(this.X);
                  this.X(i).tamara_unaries_29 = data.score{i};
               end
            end
            if this.use_logreg_11;
               data = load( TUNARIES_11, 'score' );
               for i=1:numel(this.X);
                  this.X(i).tamara_unaries_11 = data.score{i};
               end
            end
            % Shapelets
            if this.use_shapelets;
               if this.use_real_pose;
                  data = load( this.SHAPELETS_POSE, 'score' );
               else
                  data = load( this.SHAPELETS, 'score' );
               end
               for i=1:numel(this.X);
                  this.X(i).shapelets_score = data.score{i};
               end
            end
            % CPMC
            if this.use_cpmc;
               data = load( this.CPMC, 'score' );
               for i=1:numel(this.X);
                  this.X(i).cpmc_score = data.score{i};
               end
            end
            % CPMC Unaries
            if this.use_cpmc_logreg_54 || this.use_cpmc_logreg_29 || this.use_cpmc_logreg_11;
               if this.use_cpmc_logreg_54;
                  data = load( this.CPMC_UNARIES_54, 'score' );
                  for f=1:size(data.score,1);
                     for i=1:numel(this.X);
                        this.X(i).cpmc_unaries_54{f} = data.score{f,i};
                     end
                  end
               end
               if this.use_cpmc_logreg_29;
                  data = load( this.CPMC_UNARIES_29, 'score' );
                  for f=1:size(data.score,1);
                     for i=1:numel(this.X);
                        this.X(i).cpmc_unaries_29{f} = data.score{f,i};
                     end
                  end
               end
               if this.use_cpmc_logreg_11;
                  data = load( this.CPMC_UNARIES_11, 'score' );
                  for f=1:size(data.score,1);
                     for i=1:numel(this.X);
                        this.X(i).cpmc_unaries_11{f} = data.score{f,i};
                     end
                  end
               end
            end

            % Calculate similarity
            if this.use_similarity;
               if this.use_real_pose;
                  SIMILARITY = this.SIMILARITY_POSE;
               else
                  SIMILARITY = this.SIMILARITY;
               end
               data = load( SIMILARITY, 'score' );
               assert( numel(this.X) == numel(data.score) );
               for i=1:numel(this.X);
                  this.X(i).similarity_score = data.score{i};
               end
            end

            % Limb keypoints
            if this.use_limbs;
               if this.use_limbs_symm_limbs;
                  data = load( this.LIMBS, 'score' );
                  for i=1:numel(this.X);
                     this.X(i).L_score = data.score{i};
                  end
               else
                  data = load( this.JOINTS, 'score' );
                  for i=1:numel(this.X);
                     this.X(i).J_score = data.score{i};
                  end
               end

               data = load( this.LABELS, 'labels' );
               [a,b,c] = size(data.labels);
               for i=1:numel(this.X);
                  this.X(i).L_labels = reshape( data.labels(i,:,:), b, c, 1 );
               end
            end
         end

         % Load labels
         this.clothings    = sbu.Fashionista.clothings;
         % For some reason some clothes labels don't actually exist in the dataset
         % so we remove those
         if this.use_fashionista_v2; % Fashionista 0.2 has 2 more classes...
            to_remove = { 'swimwear', 'panties', 'hoodie' };
         else
            to_remove = { 'intimate', 'swimwear', 'panties', 'bodysuit', 'hoodie' };
         end
         for i=1:numel(to_remove);
            this.clothings( arrayfun( @(j) strcmp(this.clothings(j).name,to_remove{i}), 1:numel(this.clothings) ) ) = [];
         end

         % Load class distances
         % Background is missing, so that is added as a parameter.
         data    = load( sprintf(this.SIMILARITIES, this.use_disimilarity), 'wn_sim' );
         disimil = ones(size(data.wn_sim)) ./ data.wn_sim;
         disimil = disimil ./ min(disimil(:)); % Normalize to range 1-\infinity
         bg_dis  = this.loss_background * max(disimil(:));
         this.class_disimilarity = bg_dis .* ones( numel(this.clothings) );
         this.class_disimilarity(2:end,2:end) = disimil;

         % Store different settings
         this.clothings54 = this.clothings;
         this.clothings29 = this.clothings_group_fine();
         this.clothings11 = this.clothings_group_medium();

         % Group clothings to reduce the number of classes
         if strcmpi(this.class_group, 'fine')
            [clothings,disimilarity] = this.clothings_group_fine();
         elseif strcmpi(this.class_group, 'coarse')
            [clothings,disimilarity] = this.clothings_group_coarse();
         elseif strcmp(this.class_group, 'medium')
            [clothings,disimilarity] = this.clothings_group_medium();
         else
            clothings    = this.clothings;
            disimilarity = this.class_disimilarity;
         end
         this.clothings          = clothings;
         this.class_disimilarity = disimilarity;

         % Clothings maps
         this.map_54_target = this.map_clothes( this.clothings54, this.clothings );
         this.map_29_target = this.map_clothes( this.clothings29, this.clothings );
         this.map_11_target = this.map_clothes( this.clothings11, this.clothings );

         % Might as well train unaries here
         this = this.train_misc_unaries;
      end

      function [map] = map_clothes( this, clothes_in, clothes_out )
         % Special case
         if numel(clothes_in) == numel(clothes_out);
            map = eye( numel(clothes_in) );
            return
         end

         map = zeros( numel(clothes_in), numel(clothes_out) );
         for i = 1:numel(clothes_in);
            for j = 1:numel(clothes_out);
               inter = intersect( clothes_in(i).id, clothes_out(j).id );
               if (numel(inter) ~= numel(clothes_in(i).id)) || (numel(inter) ~= numel(clothes_out(j).id));
                  map(i,j) = (numel(inter) > 0);
               end
            end
         end
      end

      % Groups clothing together
      function [clothings, disimilarity] = clothings_group( this, clothes_in, disim, cgroup )
         % Build the actual group structure
         clothings   = struct( 'name', [], 'id', 0 );
         disimilarity  = zeros( numel(cgroup) );
         for i=1:numel(cgroup);
            clothings(i).name = cgroup(i).name;
            id                = [];
            for j=1:numel(cgroup(i).labels);
               id = [id,clothes_in(strcmp({clothes_in.name},cgroup(i).labels{j})).id];
            end
            clothings(i).id = id;
         end

         % Compress the disimilarity matrix
         for i=1:numel(cgroup);
            idi = ismember( [clothes_in.id], clothings(i).id );
            for j=1:numel(cgroup);
               idj = ismember( [clothes_in.id], clothings(j).id );
               dis = disim(idi,idj);
               disimilarity(i,j) = mean(dis(:));
            end
         end
      end

      % Coarse group into major groups
      function [clothings, disimilarity] = clothings_group_coarse( this, clothes_in, disim )
         if nargin < 2;
            clothes_in = this.clothings;
         end
         if nargin < 3;
            disim = this.class_disimilarity;
         end
         cgroup = struct( 'name', [], 'labels', {} );
         cgroup(1).name    = 'background';
         cgroup(1).labels  = { 'null' };
         cgroup(2).name    = 'person';
         cgroup(2).labels  = { 'hair', 'skin' };
         cgroup(3).name    = 'shoes';
         cgroup(3).labels  = { 'boots', 'clogs', 'flats', 'heels', 'loafers', 'pumps', 'sandals', 'shoes', 'socks', 'sneakers', 'wedges' };
         cgroup(4).name    = 'accessories';
         cgroup(4).labels  = { 'accessories', 'bag', 'belt', 'bracelet', 'earrings', 'glasses', 'gloves', 'hat', 'necklace', ...
                               'purse', 'ring', 'scarf', 'sunglasses', 'tie', 'wallet', 'watch' };
         cgroup(5).name    = 'upper';
         cgroup(5).labels  = { 'blazer', 'blouse', 'bodysuit', 'bra', 'cape', 'cardigan', 'coat', 'dress', 'intimate', 'jacket', 'jumper', 'romper', ...
                               'shirt', 'suit', 'sweater', 'sweatshirt', 't-shirt', 'top', 'vest' };
         cgroup(6).name    = 'lower';
         cgroup(6).labels  = { 'jeans', 'leggings', 'pants', 'shorts', 'skirt', 'stockings', 'tights' };
         % Group structure
         [clothings, disimilarity] = this.clothings_group( clothes_in, disim, cgroup );
      end

      % Fine grouping into minor groups
      function [clothings, disimilarity] = clothings_group_fine( this, clothes_in, disim )
         if nargin < 2;
            clothes_in = this.clothings;
         end
         if nargin < 3;
            disim = this.class_disimilarity;
         end
         cgroup = struct( 'name', [], 'labels', {} );
         cgroup(1).name    = 'null';
         cgroup(1).labels  = { 'null' }; % 700
         cgroup(2).name    = 'hair';
         cgroup(2).labels  = { 'hair' }; % 696
         cgroup(3).name    = 'skin';
         cgroup(3).labels  = { 'skin' }; % 700
         cgroup(4).name    = 'glasses';
         cgroup(4).labels  = {'sunglasses', 'glasses' }; % 94, 29 = 123
         cgroup(5).name    = 'accessories';
         cgroup(5).labels  = { 'accessories', 'bracelet', 'earrings', ... % 95, 72, 9
                               'necklace', 'ring', 'watch' };             % 74, 10, 32 = 292
         cgroup(6).name    = 'misc_garments';
         cgroup(6).labels  = { 'gloves', 'scarf', 'tie' }; % 12, 41, 6 = 59
         cgroup(7).name    = 'footwear';
         cgroup(7).labels  = { 'clogs', 'flats', 'heels', 'loafers', ...   % 44, 7, 70, 7
                               'pumps', 'sandals', 'sneakers', 'wedges' }; % 3, 16, 4, 55 = 206
         cgroup(8).name   = 'shoes';
         cgroup(8).labels = { 'shoes' }; % 308
         cgroup(9).name   = 'socks';
         cgroup(9).labels = { 'socks' }; % 70
         cgroup(10).name   = 'boots';
         cgroup(10).labels = { 'boots' }; % 157
         cgroup(11).name    = 'purse';
         cgroup(11).labels  = { 'purse', 'wallet' }; % 65, 9 = 74
         cgroup(12).name   = 'bag';
         cgroup(12).labels = { 'bag' }; % 255
         cgroup(13).name    = 'hat';
         cgroup(13).labels  = { 'hat' }; % 104
         cgroup(14).name   = 'dress';
         cgroup(14).labels = { 'dress', 'romper' }; % 192, 3 = 195
         cgroup(15).name   = 'upper_body';
         cgroup(15).labels = { 'cardigan', 'coat', 'suit', 'bodysuit' }; % 72, 68, 3 = 197 (bodysuit??)
         cgroup(16).name   = 'shirt';
         cgroup(16).labels = { 'shirt', 't-shirt' }; % 86, 72 = 158
         cgroup(17).name   = 'sweater';
         cgroup(17).labels = { 'jumper', 'sweater', 'sweatshirt' }; % 31, 58, 2 = 91
         cgroup(18).name   = 'top';
         cgroup(18).labels = { 'bra', 'top' }; % 3, 147 = 150
         cgroup(19).name   = 'jacket';
         cgroup(19).labels = { 'blazer', 'jacket' }; % 80, 83 = 163
         cgroup(20).name   = 'vest';
         cgroup(20).labels = { 'vest', 'cape' }; % 3, 51 = 54
         cgroup(21).name   = 'blouse';
         cgroup(21).labels = { 'blouse' }; % 120
         cgroup(22).name   = 'jeans';
         cgroup(22).labels = { 'jeans' }; % 86
         cgroup(23).name   = 'leggings';
         cgroup(23).labels = { 'leggings' }; % 63
         cgroup(24).name   = 'pants';
         cgroup(24).labels = { 'pants' }; % 94
         cgroup(25).name   = 'shorts';
         cgroup(25).labels = { 'shorts', 'intimate' }; % 113 (intimate??)
         cgroup(26).name   = 'skirt';
         cgroup(26).labels =  { 'skirt' }; % 128
         cgroup(27).name   = 'stockings';
         cgroup(27).labels = { 'stockings' }; %63
         cgroup(28).name   = 'tights';
         cgroup(28).labels = { 'tights' }; % 87
         cgroup(29).name   = 'belt';
         cgroup(29).labels = { 'belt' }; % 185
         % Group structure
         [clothings, disimilarity] = this.clothings_group( clothes_in, disim, cgroup );
      end

      % Fine grouping into minor groups
      function [clothings, disimilarity] = clothings_group_medium( this, clothes_in, disim )
         if nargin < 2;
            clothes_in = this.clothings;
         end
         if nargin < 3;
            disim = this.class_disimilarity;
         end
         cgroup = struct( 'name', [], 'labels', {} );
         cgroup(1).name    = 'null';
         cgroup(1).labels  = { 'null' }; % 700
         cgroup(2).name    = 'hair';
         cgroup(2).labels  = { 'hair' }; % 696
         cgroup(3).name    = 'skin';
         cgroup(3).labels  = { 'skin' }; % 700
         cgroup(4).name    = 'glasses';
         cgroup(4).labels  = {'sunglasses', 'glasses' }; % 94, 29 = 123
         cgroup(5).name    = 'accessories';
         cgroup(5).labels  = { 'accessories', 'bracelet', 'earrings', ... % 95, 72, 9
                               'gloves', 'scarf', 'tie', ... % 12, 41, 6
                               'necklace', 'ring', 'watch' };             % 74, 10, 32 = 292
         cgroup(6).name    = 'feeties';
         cgroup(6).labels  = { 'clogs', 'flats', 'heels', 'loafers', ...   % 44, 7, 70, 7
                               'pumps', 'sandals', 'sneakers', 'wedges', ...% 3, 16, 4, 55
                               'shoes', 'socks', 'boots' };               % 308, 70, 157
         cgroup(7).name    = 'bag';
         cgroup(7).labels  = { 'bag', 'purse', 'wallet' }; % 255, 65, 9
         cgroup(8).name    = 'hat';
         cgroup(8).labels  = { 'hat' }; % 104
         cgroup(9).name   = 'torso';
         cgroup(9).labels = { 'bodysuit', 'cardigan', 'coat', 'suit', ... % 72, 68, 3 (bodysuit??)
                               'shirt', 't-shirt', ... 86, 72
                               'jumper', 'sweater', 'sweatshirt', ... % 31, 58, 2 = 91
                               'bra', 'top', ... % 3, 147
                               'blazer', 'jacket', ... % 80, 83
                               'vest', 'cape', ... % 3, 51 
                               'blouse', ... % 120
                               'dress', 'romper' }; ... % 192, 3 
         cgroup(10).name   = 'legs';
         cgroup(10).labels = { 'jeans', 'leggings', 'pants', 'shorts', 'skirt', 'stockings', 'tights', 'intimate' }; % 86, 63, 94, 113, 128, 63, 87 (intimate??)
         cgroup(11).name   = 'belt';
         cgroup(11).labels = { 'belt' }; % 185
         % Group structure
         [clothings, disimilarity] = this.clothings_group( clothes_in, disim, cgroup );
      end


      % Trains all the features needed in the proper order
      function this = train_all(this)
         this = this.train_misc_unaries;
         this = this.train_MRF;
      end
      
      function x = path(this, file_name)
         %GET.PATH storage file path
         if nargin < 2, file_name = ''; end
         x = fullfile( this.TMP, this.name, file_name );
         if ~exist(fileparts(x),'dir'), mkdir(fileparts(x)); end
      end

      function eval_details( this, R )
         fprintf( '%16s   Accu.  Prec.  Reca.  VOC   Occurences\n', '' );
         fprintf( '%16s   %4.2f   %4.2f   %4.2f   %4.2f      \n', 'All', ...
                  R.acc, R.avr_pre, R.avr_rec, R.avr_voc );
         [foo,order] = sort(this.class_occurences, 2, 'descend');
         for i=order;
            fprintf( '%16s          %4.2f   %4.2f   %4.2f   %.1e\n', ...
                  this.clothings(i).name, ...
                  R.pre(i), R.rec(i), R.voc(i), this.class_occurences(i) );
         end
      end

      % Evaluates the segmentation for a part on the annotated ground truth
      function [R] = eval_segmentation( this, ids, Lpred )
         if nargin < 2;
            ids = this.testids;
         end

         Ctrue = cell( numel(ids), 1 );
         Cpred = cell( numel(ids), 1 );
         for i=1:numel(ids);
            id = ids(i);
            ts = tic;
            fprintf( '\r%03d / %03d... ', i, numel(ids) );
            gt       = load( sprintf(this.GTMAPS, id), 'gtmap' );
            Gmap     = gt.gtmap(:);
            Ctrue{i} = nan(size(Gmap));
            for j=1:numel(this.clothings);
               Ctrue{i}(ismember(Gmap, this.clothings(j).id)) = j;
            end
            Ctrue{i} = uint8(Ctrue{i});
            seg      = this.photos_truth(id).segmentation;
            pred     = seg.map;
            for j=1:seg.length;
               pred( seg.map==j ) = Lpred{i}(j);
            end
            Cpred{i} = uint8(pred(:));
         end
         fprintf( '\n' );
         % Confusion matrix and order
         [C,O] = confusionmat( cat(1,Ctrue{:}), cat(1,Cpred{:}), 'Order', (1:numel(this.clothings))' );
         R           = struct;
         R.confusion = C;
         R.order     = O;
         % Per class
         tp    = diag(C);
         fp    = sum(C,1)'; % includes tp
         fn    = sum(C,2); % includes tp
         R.acc = sum(tp) / sum(C(:));
         R.pre = tp ./ fp;
         R.rec = tp ./ fn;
         R.f1  = 2*tp ./ (fn+fp);
         R.voc = tp ./ (fp+fn-tp);
         % Means
         R.avr_pre = mean(R.pre(~isnan(R.pre)));
         R.avr_rec = mean(R.rec(~isnan(R.rec)));
         R.avr_f1  = mean(R.f1( ~isnan(R.f1)));
         R.avr_voc = mean(R.voc(~isnan(R.voc)));
      end

      % Evaluates oracle performance (should be 100% when using ground truth superpixels)
      function [R] = test_oracle( this, ids )
         if nargin < 2;
            ids = 1:numel(this.photos_truth);
         end
         L = cell( numel(ids), 1 );
         for i=1:numel(ids);
            data = load( sprintf(this.SEGMENTS, ids(i)), 'labels' );
            L{i} = nan(size(data.labels));
            for j=1:numel(this.clothings);
               L{i}(ismember(data.labels, this.clothings(j).id)) = j;
            end
         end
         R = this.eval_segmentation( ids, L );
      end

      % Exports the MRF so it can be run with mdSP
      function MRF_export( this, path )
         real_pose = this.use_real_pose;

         if ~exist( path, 'dir' ); mkdir( path ); end
         pathtrain = [path,'/Train/'];
         pathtest  = [path,'/Test/'];
         if ~exist( pathtrain, 'dir' ); mkdir( pathtrain ); end
         if ~exist( pathtest, 'dir' );  mkdir( pathtest );  end
         % Train
         if real_pose;
            %this.use_real_pose = 0;
            %this = this.features_load;
         end
         [S,w,wlabels] = this.build_MRF( this.trainids );
         this.wlabels = wlabels'; % Save labels
         WriteOutput( pathtrain, S );
         % Test
         if real_pose;
            %this.use_real_pose = 1;
            %this = this.features_load;
         end
         S = this.build_MRF( this.testids );
         WriteOutput( pathtest,  S );
         % Model parameters
         model = this.struct;
         model.X = [];
         model.photos_truth = [];
         save( [path,'/model.mat'], 'model' );
      end

      % Visualizes weight values
      function MRF_weights_view( this )
         arrayfun( @(i) fprintf('%f %s\n',this.weights(i),this.wlabels{i}), 1:numel(this.weights) );
      end

      % Trains and tests the MRF on the same samples
      function [R] = traintest_MRF( this, ids )
         if nargin < 2;
            ids = this.trainids;
         end
         [S, w] = this.build_MRF( ids );
         ts = tic;
         [w, Pred] = structuredPrediction( S, ...
               this.Params, 1, w );
         fprintf( 'Trained and tested MRF in %.1f seconds\n', toc( ts ) );

         NUM_CLASSES = numel(this.clothings);
         Lpred    = cell( numel(Pred), 1 );
         for i=1:numel(Pred);
            regions  = cat( 1, Pred{i}.Regions{:} );
            beliefs  = cat( 2, regions(arrayfun( @(j) size(regions(j).Beliefs,1), 1:numel(regions) ) == NUM_CLASSES ).Beliefs );
            [foo,b]  = max( beliefs );
            Lpred{i} = b';
         end
         R = this.eval_segmentation( ids, Lpred );
      end

      function test_MRF_visualize( this, path, Pred, ids )
         if nargin < 3;
            Pred = this.test_MRF();
         end
         if nargin < 4;
            ids = this.testids;
         end
         mkdir( path );
         NUM_CLASSES = numel(this.clothings);
         for i=1:numel(Pred);
            fprintf( '\r%03d / %03d', i, numel(Pred ) );
            p        = this.photos_truth(ids(i)).copy;
            seg      = p.segmentation;
            nsegments = seg.length;
            % Get labels
            regions  = cat( 1, Pred{i}.Regions{1:nsegments} );
            beliefs  = cat( 2, regions(arrayfun( @(j) size(regions(j).Beliefs,1), 1:numel(regions) ) == NUM_CLASSES ).Beliefs );
            [foo,b]  = max( beliefs );
            Lpred    = b';
            pred_label = ones(numel(this.clothings),1);
            pred_label(unique(b)) = 2;

            % Visualize
            FNAME = sprintf( '%s/%06d', path, ids(i) );
            if ~exist( [FNAME,'_gt.png'], 'file' );
               %seg.show(  'image', p.image, 'labels', this.clothings );
               seg.show(  'labels', this.clothings );
               colorbar off;
               export_fig( [FNAME,'_gt.png'], '-png', '-transparent' );
            end
            if ~exist( [FNAME,'_pred.png'], 'file' );
               segL  = seg.copy();
               segL.labels = arrayfun( @(i) this.clothings(i).id(1), Lpred );
               %segL.show( 'image', p.image, 'labels', this.clothings );
               segL.show( 'labels', this.clothings );
               colorbar off;
               export_fig( [FNAME,'_pred.png'], '-png', '-transparent' );
            end
         end
      end

      function [R] = test_MRF_segmentation( this, Pred, ids )
         if nargin < 2;
            Pred = this.test_MRF();
         end
         if nargin < 3;
            ids = this.testids;
         end
         NUM_CLASSES = numel(this.clothings);
         Lpred       = cell( numel(Pred), 1 );
         for i=1:numel(Pred);
            nsegments = this.photos_truth(ids(i)).segmentation.length;
            regions  = cat( 1, Pred{i}.Regions{1:nsegments} );
            beliefs  = cat( 2, regions(arrayfun( @(j) size(regions(j).Beliefs,1), 1:numel(regions) ) == NUM_CLASSES ).Beliefs );
            [foo,b]  = max( beliefs );
            Lpred{i} = b';
         end
         R = this.eval_segmentation( ids, Lpred );
      end

      % Tests the MRF
      function [Pred] = test_MRF( this, ids )
         if nargin < 2;
            ids = this.testids;
         end
         S     = this.build_MRF( ids );
         ts = tic;
         [w, Pred] = structuredPrediction( S, ...
               this.Params, 2, this.weights );
         fprintf( 'Tested MRF in %.1f seconds\n', toc( ts ) );
      end

      % Trains the MRF
      function this = train_MRF( this, ids )
         if nargin < 2;
            ids = this.trainids;
         end
         [S,w_init,wlabels] = this.build_MRF( ids );
         this.wlabels = wlabels';
         if numel(this.weights)>0; % Use as a reference if it's set already
            n           = min(numel(w_init),numel(this.weights));
            w_init(1:n) = this.weights(1:n);
         end
         ts = tic;
         fprintf( 'Weights to learn = %d\n', numel(w_init) );
         [this.weights, P] = structuredPrediction( ...
               S, this.Params, 0, w_init );
         fprintf( 'Trained MRF in %.1f seconds\n', toc( ts ) );
      end

      % Builds the MRF from photos and X
      function [S,w,wlabels] = build_MRF( this, ids )
         % Display Info
         fprintf( 'BUILDING MRF OUTPUT %d CLASSES (REAL POSE=%d)...\n', numel(this.clothings), this.use_real_pose );
         fprintf( 'UNARIES:\n' );
         if this.use_fullbias;
            fprintf( '   fullbias\n' );
         else
            fprintf( '   bgbias\n' );
         end
         if this.use_logreg_54 || this.use_logreg_29 || this.use_logreg_11;
            fprintf( '   logreg:     ' );
            if this.use_logreg_54; fprintf( '  54' ); end
            if this.use_logreg_29; fprintf( '  29' ); end
            if this.use_logreg_11; fprintf( '  11' ); end
            fprintf( '\n' );
         end
         if this.use_cpmc_logreg_54 || this.use_cpmc_logreg_29 || this.use_cpmc_logreg_11;
            fprintf( '   cpmc_logreg:' );
            if this.use_cpmc_logreg_54; fprintf( '  54' ); end
            if this.use_cpmc_logreg_29; fprintf( '  29' ); end
            if this.use_cpmc_logreg_11; fprintf( '  11' ); end
            fprintf( '\n' );
         end
         if this.use_cpmc;      fprintf( '   cpmc\n' ); end
         if this.use_shapelets; fprintf( '   shapelets\n' ); end
         fprintf( 'HIGHER ORDER\n' );
         if this.use_similarity;fprintf( '   similarity\n' ); end
         if this.use_limbs;fprintf(      '   limbs\n' ); end

         photos = this.photos_truth( ids );
         X      = this.X( ids );
         
         O              = struct;
         % Global stuff
         O.NUM_IMAGES   = numel(photos); % Test on a few
         O.NUM_CLASSES  = numel(this.clothings);
         O.NUM_CLOTHES  = O.NUM_CLASSES;

         % Start loading the structure
         for i=1:O.NUM_IMAGES;
            numvars  = 0; % Number of variables counter
            ts       = tic;
            fprintf( '\rInitializing Image %03d / %03d...', i, O.NUM_IMAGES );
            seg         = photos(i).segmentation;
            O.nsegments = numel(seg.labels);

            % Each segment gets it's own random variable with as many states
            % as we have segmantic classes
            S{i}.VariableCardinalities = O.NUM_CLASSES.*ones(O.nsegments,1); % Segmentation
            if this.use_limbs;
               if this.use_limbs_symm_limbs;
                  limbs = this.symm_limbs;
               else
                  limbs = this.symm_joints;
               end
               n = numel( S{i}.VariableCardinalities );
               if this.use_limbs_merge;
                  S{i}.VariableCardinalities = [S{i}.VariableCardinalities;...
                        O.NUM_CLOTHES.*ones(size(limbs,1),1)]; % Zs
                  IdsL = reshape( n+(1:(size(limbs,1))), ...
                                  size(limbs,1), 1 );
               else
                  S{i}.VariableCardinalities = [S{i}.VariableCardinalities;...
                        O.NUM_CLOTHES.*ones(size(limbs,1)*2,1)]; % Zs
                  IdsL = reshape( n+(1:(size(limbs,1)*2)), ...
                                  size(limbs,1), 2 );
               end

               % Load coocurrence
               data = load( this.LABELS );
               L_cooc = this.L_coocurrence( data.labels );
            end
            S{i}.Observation  = nan( size(S{i}.VariableCardinalities) );
            S{i}.MachineIDs   = (mod( i, this.NUM_MACHINES)+1) .* ones( size(S{i}.VariableCardinalities) );

            % Build segments
            wlabels = {};
            [S{i}, numvars, IdsSeg, wlabels] = this.build_MRF_segments( S{i}, O, seg, numvars, X(i), wlabels );
            if this.use_similarity;
               [S{i}, numvars, wlabels] = this.build_MRF_similarity( S{i}, O, seg, numvars, X(i), IdsSeg, wlabels );
            end
            if this.use_limbs;
               if this.use_limbs_merge;
                  [S{i}, numvars, wlabels] = this.build_MRF_seg_symmetry_merge( S{i}, O, seg, numvars, X(i), IdsSeg, IdsL, wlabels, L_cooc );
               else
                  [S{i}, numvars, wlabels] = this.build_MRF_seg_symmetry( S{i}, O, seg, numvars, X(i), IdsSeg, IdsL, wlabels, L_cooc );
               end
            end
            fprintf( '   %.1f seconds!   ', toc(ts) );
         end
         fprintf( '\n' );
         w = ones( numvars, 1 ) + 1e-3 .* rand( numvars, 1 ); % Break possibly symmetry issues with some randomness
      end

      function [this, R, scores] = compute_cpmc_logreg_unaries_clothings( this, clothings, trainids, testids )
         if nargin < 3;
            trainids = this.trainids;
            testids  = this.testids;
         end

         featpath = this.CPMCFEATS;
         featlist =  {'CPMC_segms_150_sp_approx_LBP_f_pca_2500_noncent', ...
                      'CPMC_segms_150_sp_approx_SIFT_GRAY_f_g_pca_5000_noncent', ...
                      'CPMC_segms_150_sp_approx_SIFT_GRAY_mask_pca_5000_noncent'};

         for i=1:numel(this.photos_truth);
            this.X(i).cpmc_unaries = cell( numel(featlist), 1 );
         end

         scores   = cell( numel(featlist), numel(this.photos_truth) );
         R        = cell( numel(featlist), 1 );
         for f=1:numel(featlist);
            % TRAINING
            L = cat(1,this.X(trainids).L);
            F = [];
            for i = trainids;
               data  = load( sprintf( '%s%s/%06d.mat', featpath, featlist{f}, i ), 'D' );
               F     = [F; data.D'];
            end
            F = double(F);
            logreg = cell( numel(clothings), 1 );
            for i=1:numel(clothings);
               fprintf( '\rTraining Logistic Regression %02d / %02d for %s    ', ...
                     i, numel(clothings), featlist{f} );
               logreg{i} = liblinear.SVM( ismember(L,clothings(i).id), F, ...
                     'prob', true, 'type', 0, 'bias', 1, 'autoweight', true, 'quiet', true );
            end
            fprintf( '\n' );

            % COMPUTE AS FEATURE
            for i = 1:numel(this.photos_truth);
               L     = cat(1,this.X(i).L);
               data  = load( sprintf( '%s%s/%06d.mat', featpath, featlist{f}, i ), 'D' );
               F     = double( data.D' );
               rescell  = arrayfun( @(j) max([zeros(numel(L),1) ...
                     logreg{j}.predict( L, F ).val(:,logreg{j}.model.Label==1)], ...
                  [],2), 1:numel(clothings), 'UniformOutput', false );
               this.X(i).cpmc_unaries{f} = cat( 2, rescell{:} );
               scores{f,i} = this.X(i).cpmc_unaries{f};
            end

            % EVALUATION
            results = [];
            for i = testids;
               results  = [results;this.X(i).cpmc_unaries{f}];
            end
            [foo,P]  = max( results, [], 2 );

            I  = arrayfun( @(j) j*ones( size(this.X(j).L) ), testids, 'UniformOutput', false );
            I  = cat(1,I{:});
            L  = cell( numel(testids), 1 );
            for i=1:numel(testids);
               L{i} = P(I==testids(i));
            end
         end
      end

      function this = compute_cpmc_logreg_unaries( this, path )
         TUNARIES_54 = this.CPMC_UNARIES_54_FILE;
         TUNARIES_29 = this.CPMC_UNARIES_29_FILE;
         TUNARIES_11 = this.CPMC_UNARIES_11_FILE;
         if ~exist( [path TUNARIES_54], 'file' );
            [a,b,score] = this.compute_cpmc_logreg_unaries_clothings( this.clothings54 );
            save( [path TUNARIES_54], 'score' );
         end
         if ~exist( [path TUNARIES_29], 'file' );
            [a,b,score] = this.compute_cpmc_logreg_unaries_clothings( this.clothings29 );
            save( [path TUNARIES_29], 'score' );
         end
         if ~exist( [path TUNARIES_11], 'file' );
            [a,b,score] = this.compute_cpmc_logreg_unaries_clothings( this.clothings11 );
            save( [path TUNARIES_11], 'score' );
         end
      end

      function logreg = this.logreg( this )
         if strcmpi(this.class_group, 'fine')
            logreg = this.logreg29;
         elseif strcmp(this.class_group, 'medium')
            logreg = this.logreg11;
         else
            logreg = this.logreg54;
         end
      end

      function this = compute_logreg_unaries( this, path )
         if this.use_real_pose;
            TUNARIES_54 = this.TAMARA_UNARIES_POSE_54_FILE;
            TUNARIES_29 = this.TAMARA_UNARIES_POSE_29_FILE;
            TUNARIES_11 = this.TAMARA_UNARIES_POSE_11_FILE;
         else
            TUNARIES_54 = this.TAMARA_UNARIES_54_FILE;
            TUNARIES_29 = this.TAMARA_UNARIES_29_FILE;
            TUNARIES_11 = this.TAMARA_UNARIES_11_FILE;
         end
         if ~exist( [path TUNARIES_54], 'file' );
            score = this.compute_logreg_unaries_clothings( this.clothings54 );
            save( [path TUNARIES_54], 'score' );
         end
         if ~exist( [path TUNARIES_29], 'file' );
            score = this.compute_logreg_unaries_clothings( this.clothings29 );
            save( [path TUNARIES_29], 'score' );
         end
         if ~exist( [path TUNARIES_11], 'file' );
            score = this.compute_logreg_unaries_clothings( this.clothings11 );
            save( [path TUNARIES_11], 'score' );
         end
      end

      function score = compute_logreg_unaries_clothings( this, clothings )
         % Load full data
         F      = load( this.FEATURES, 'Xtruth' );
         X      = F.Xtruth;

         % Set up features
         L = cat(1,X(this.trainids).L);
         F = cat(1,X(this.trainids).F);

         % Train features on liblinear
         logreg = cell( numel(clothings), 1 );
         ts = tic;
         for i=1:numel(clothings);
            fprintf( '\rTraining Logistic Regression %02d / %02d    ', ...
                  i, numel(clothings) );
            logreg{i} = liblinear.SVM( ismember(L,clothings(i).id), F, ...
                  'prob', true, 'type', 0, 'bias', 1, 'autoweight', true, 'quiet', true );
         end
         fprintf( '\nTrained Logistic Regression in %.1f seconds.\n', toc(ts) );

         % Evaluate
         if this.use_real_pose;
            F      = load( this.FEATURES_POSE, 'Xtruth' );
         end

         score{i} = cell( numel(this.photos_truth) );
         for i = 1:numel(this.photos_truth);
            L     = cat(1,this.X(i).L);
            seg = this.photos_truth(i).segmentation;
            F     = X(i).F;
            rescell     = arrayfun( @(j) max([zeros(seg.length,1) ...
                     logreg{j}.predict( L, F ).val(:,logreg{j}.model.Label==1)], ...
                  [],2), 1:numel(clothings), 'UniformOutput', false );
            score{i} = cat( 2, rescell{:} ); % Results of running the logistic regressions
         end
      end

      function R = test_logreg_unaries( this, ids )
         if nargin < 2;
            ids = this.testids;
         end
         logreg = this.logreg();
         % Test evaluation
         I  = arrayfun( @(j) j*ones( size(this.X(j).L) ), ids, 'UniformOutput', false );
         I  = cat(1,I{:});
         L  = cat(1,this.X(ids).L);
         F  = cat(1,this.X(ids).F);
         ts = tic;
         rescell  = arrayfun( @(j) max([zeros(numel(L),1) logreg{j}.predict( L, F ).val(:,logreg{j}.model.Label==1)],[],2), 1:numel(this.clothings), 'UniformOutput', false );
         fprintf( 'Tested Logistic Regression in %.1f seconds.\n', toc(ts) );
         results  = cat( 2, rescell{:} );
         [foo,P]  = max( results, [], 2 );
         % Organize output
         L = cell( numel(ids), 1 );
         for i=1:numel(ids);
            L{i} = P(I==ids(i));
         end
         R = this.eval_segmentation( ids, L );
      end

      function this = train_misc_unaries( this, ids )
         if nargin < 2;
            ids      = this.trainids;
         end

         % Probabilities
         O = zeros(1,numel(this.clothings));
         P = zeros(1,numel(this.clothings));
         for i=ids;
            P = P + arrayfun( @(j) (sum(ismember(this.X(i).L, this.clothings(j).id))>0), 1:numel(this.clothings) );
            O = O + arrayfun( @(j) sum(this.X(i).area(find(ismember(this.X(i).L,this.clothings(j).id)))), 1:numel(this.clothings) );
         end
         this.class_occurences   = O ./ numel(this.trainids);
      end

      function C = L_coocurrence( this, labels, ids )
         if nargin < 3;
            ids      = this.trainids;
         end

         % Calculate correlation matrix
         NUM_CLASSES = numel(this. clothings);
         C = zeros( NUM_CLASSES, NUM_CLASSES, size(this.symm_limbs_tree,1), 2 );
         L = labels( ids, :, : );
         for p = 1:size(this.symm_limbs_tree,1);
            i1 = this.symm_limbs_tree(p,1);
            i2 = this.symm_limbs_tree(p,2);
            for i=1:numel(this.clothings);
               C(:,i,p,1) = arrayfun( @(l) sum( L( L(:,i1,1)==i, i2, 1 ) == l ), 1:NUM_CLASSES )';
               C(:,i,p,2) = arrayfun( @(l) sum( L( L(:,i1,2)==i, i2, 2 ) == l ), 1:NUM_CLASSES )';
            end
         end

         % Normalize for samples
         C = C ./ numel(ids);
      end

      % Computes shapelet features
      function [features,shapelets] = shapelets_compute_features( this )
         [shapelets, Z] = this.shapelets_compute;
         features = this.shapelets_evaluate( this.photos_truth, shapelets );
      end

      % Gets the amount of shapelet on each segment
      function [scores] = shapelets_train_size( this )
         NUM_CLASSES = numel(this.clothings);
         sha_size    = 100;
         scalelist   = 0.1:0.1:3.0;
         %scalelist   = 0.5:0.1:2.0;
         ids         = this.trainids;

         scores = zeros( numel(scalelist), NUM_CLASSES, max(ids), 3 );
         for s = 1:numel(scalelist);
            sha_scale = scalelist(s);

            cachefile = sprintf( 'cache/size%d_scale%.1f.mat', sha_size, sha_scale );
            if exist( cachefile, 'file' );
               fprintf('Loading Shapelets for scale = %.1f\n', sha_scale );
               data           = load( cachefile, 'score_data' );
               scores(s,:,:,:) = data.score_data;
               continue;
            end

            shapelets = this.shapelets_compute( ids, sha_size, sha_scale );
            for i = ids;
               fprintf( '\rTraining Shapelets %03d / %03d for scale = %.1f   ', find(i==ids), numel(ids), sha_scale );
               seg   = this.photos_truth(i).segmentation;
               smap  = seg.label_map;
               UCI   = this.photos_truth(i).to_uci;
               S     = this.shapelets_map( UCI, shapelets );
               S     = S{1};

               for j=1:NUM_CLASSES;
                  M  = ismember( smap, this.clothings(j).id);
                  TP = S(:,:,j) .*  M;
                  FP = S(:,:,j)  - TP;
                  FN =    M      - TP;
                  scores( s, j, i, 1 ) = sum(TP(:));
                  scores( s, j, i, 2 ) = sum(FP(:));
                  scores( s, j, i, 3 ) = sum(FN(:));
               end
            end
            fprintf('\n');
            score_data = scores( s,:,:,: );
            save( cachefile, 'score_data' );
         end

         ss = mean( scores(:,:,ids,1) ./ sum(scores(:,:,ids,:),4), 3 );
         for j=1:NUM_CLASSES;
            [foo,i] = max( ss(:,j) );
            fprintf( 'Shapelets %16s, best result is %.1f\n', this.clothings(j).name, scalelist(i) );
         end
      end

      % Computes shapelets per clothings
      function [shapelets, Z] = shapelets_compute( this, ids, sha_size, sha_scale )
         if nargin < 2;
            ids = this.trainids;
         end
         if nargin < 3;
            sha_size  = 100;
            sha_scale = this.shape_size;
         end

         % Set up stuff
         UCI            = this.photos_truth.to_uci;
         NUM_PARTS      = numel(this.pa);
         NUM_CLASSES    = numel(this.clothings);
         SHAPELET_SIZE  = sha_size;
         %SHAPELET_SCALE = sha_scale;
         if numel(sha_scale==1);
            sha_scale = repmat( sha_scale, 1, NUM_CLASSES );
         end
       
         % Compute the shapelets
         shapelets = zeros( SHAPELET_SIZE, SHAPELET_SIZE, NUM_CLASSES, NUM_PARTS );
         Z         = zeros( NUM_CLASSES, 1 ); % Appearences (to normalize)
         for i=ids;
            fprintf( '\rShapelets %03d / %03d   ', find(ids==i), numel(ids) );
            for j=1:NUM_CLASSES;
               SHAPELET_SCALE = this.shape_size(j);

               % Ignore if not in the image
               if ~ismember(this.photos_truth(i).segmentation.labels, this.clothings(j).id); continue; end;
               Z(j) = Z(j)+1;

               box = pointtobox( UCI(i), this.pa );
               I   = ismember( UCI(i).labels, this.clothings(j).id );
               for k=1:NUM_PARTS;
                  x1 = box.x1(k); x2 = box.x2(k);
                  y1 = box.y1(k); y2 = box.y2(k);
                  dx = x2 - x1;   dy = y2 - y1;
                  ox = (SHAPELET_SCALE-1)*dx/2;
                  oy = (SHAPELET_SCALE-1)*dy/2;
                  x1 = x1-ox;     x2 = x2+ox;
                  y1 = y1-oy;     y2 = y2+oy;
                  IC = imcrop( I, [x1, y1, x2-x1, y2-y1] );
                  IR = imresize( IC, [SHAPELET_SIZE, SHAPELET_SIZE] );
                  shapelets(:,:,j,k) = shapelets(:,:,j,k) + IR;
               end
            end
         end
         for j=1:NUM_CLASSES;
            shapelets(:,:,j,:) = shapelets(:,:,j,:) ./ Z(j);
         end
         %shapelets = shapelets ./ numel(ids); % Normalize
         fprintf('\n');
      end

      % Gets the amount of shapelet on each segment
      function [shapelets_score] = shapelets_evaluate( this, photos, shapelets )
         NUM_CLASSES       = size(shapelets,3);
         shapelets_score   = cell( numel(photos), 1 );
         for i=1:numel(photos);
            fprintf( '\rEvaluating Shapelets %03d / %03d   ', i, numel(photos) );
            seg   = photos(i).segmentation;
            nsegments = seg.length;
            UCI   = photos(i).to_uci;
            I     = this.shapelets_map( UCI, shapelets );
            score = zeros( nsegments, NUM_CLASSES );
            for j=1:nsegments;
               for k=1:NUM_CLASSES;
                  Ik          = I{1}(:,:,k);
                  score(j,k)  = sum(Ik(seg.map==j));
               end
               score(j,:) = score(j,:) ./ seg.area(j); % Normalize
            end
            shapelets_score{i} = score;
         end
         fprintf('\n');
      end

      function [shapelets_score] = shapelets_evaluate_gt( this, photos )
         NUM_CLASSES       = numel(this.clothings);
         shapelets_score   = cell( numel(photos), 1 );
         for i=1:numel(photos);
            fprintf( '\rEvaluating Shapelets %03d / %03d   ', i, numel(photos) );
            seg   = photos(i).segmentation;
            nsegments = seg.length;
            score = zeros( nsegments, NUM_CLASSES );
            gt    = load( sprintf(this.GTMAPS, i), 'gtmap' );
            for j=1:nsegments;
               for k=1:NUM_CLASSES;
                  score(j,k) = sum(ismember(gt.gtmap(seg.map==j), ...
                                            this.clothings(k).id));
               end
               score(j,:) = score(j,:) ./ seg.area(j); % Normalize
            end
            shapelets_score{i} = score;
         end
         fprintf('\n');
      end

      % Evaluates a shapelet on a set of images.
      % Produces a set of heat maps corresponding to the shapelet output.
      function [I] = shapelets_map( this, UCI, shapelets, sha_size, sha_scale )
         if nargin < 4;
            sha_size = 100;
            sha_scale = this.shape_size;
         end
         %SHAPELET_SCALE = sha_scale;
         NUM_PARTS   = size(shapelets,4);
         NUM_CLASSES = size(shapelets,3);
         SHAPELET_SIZE = sha_size;
         if numel(sha_scale==1);
            sha_scale = repmat( sha_scale, 1, NUM_CLASSES );
         end
         [c,h]       = size(UCI(1).labels);
         I           = cell( numel(UCI), 1 );
         for i=1:numel(UCI);
            Iuci  = zeros( c, h, NUM_CLASSES, NUM_PARTS );
            %Iuci(:,:,1,:) = 1;
            bgmask = zeros( c, h, NUM_PARTS );
            box   = pointtobox( UCI(i), this.pa );
            for k=1:NUM_PARTS;
               x1 = box.x1(k); x2 = box.x2(k);
               y1 = box.y1(k); y2 = box.y2(k);
               dx = x2 - x1;   dy = y2 - y1;

               %IO = zeros( c, h, NUM_CLASSES );
               for j=1:NUM_CLASSES;
                  SHAPELET_SCALE = this.shape_size(j);

                  ox = (SHAPELET_SCALE-1)*dx/2;
                  oy = (SHAPELET_SCALE-1)*dy/2;
                  x1 = x1-ox;     x2 = x2+ox;
                  y1 = y1-oy;     y2 = y2+oy;
                  fx = floor(x1); fy = floor(y1);
                  cx = ceil(x2);  cy = ceil(y2);
                  dxf = (cx-fx); dxp = (x2-x1);
                  dyf = (cy-fy); dyp = (y2-y1);
                  xr = fx:cx;     yr = fy:cy;
                  ss2 = (SHAPELET_SIZE+1)/2;
                  [xs,ys] = meshgrid( (1:SHAPELET_SIZE), (1:SHAPELET_SIZE) );
                  [xq,yq] = meshgrid( (xr-mean(xr))/dxf*dxp/dxf*SHAPELET_SIZE+ss2, ...
                                      (yr-mean(yr))/dyf*dxp/dxf*SHAPELET_SIZE+ss2 );

                  z  = interp2( xs, ys, shapelets(:,:,j,k), xq, yq, 'linear', nan );
                  px = (xr>0) & (xr<=h);
                  py = (yr>0) & (yr<=c);
                  Iuci(yr(py),xr(px),j,k) = z(yr(py)-fy+1,xr(px)-fx+1);
                  if j==1;
                     bgmask(yr(py),xr(px),k) = ~isnan(z(yr(py)-fy+1,xr(px)-fx+1));
                  end
               end
            end
            bg   = (sum(bgmask,3)==0); 
            Iuci(isnan(Iuci)) = 0;
            I{i} = sum( Iuci, 4 ) ./ sum( Iuci>0, 4 ); % Compress parts together with max
            I{i}(sum(Iuci>0,4)==0) = 0; % Remove spurious nans
            Ibg  = I{i}(:,:,1);
            Ibg(bg) = 1;
            I{i}(:,:,1) = Ibg;
         end
      end

      function symmetry_visualize( this )
         w1 = this.symm_extension;
         w2 = this.symm_width;
         limbs = this.symm_limbs;
         for p=1:numel(this.photos_truth);
            score = this.X(p).L_score;
            limbs = this.symm_limbs;
            seg = this.photos_truth(p).segmentation;

            %imshow(I ./ max(I(:)))
            figure(1);
            clf;
            seg.show_superpixels;
            hold on;
            col = distinguishable_colors( size(limbs,1) );
            u = zeros(seg.length,1);
            v = zeros(seg.length,1);
            for s = 1:seg.length;
               [I,J] = find( seg.map==s );
               u(s) = mean(I);
               v(s) = mean(J);
            end
            n = 0;

            for i = 1:size(limbs,1);
               inrange1 = find( score(:,i,1) > this.symm_thresh )';
               inrange2 = find( score(:,i,2) > this.symm_thresh )';
               if numel(inrange1)==0 || numel(inrange2)==0; continue; end
               for i1=inrange1;
                  for i2=inrange2;
                     if i1==i2; continue; end
                     n = n+1;
                     ss = score(i1,i,1) * score(i2,i,2);
                     plot( v([i1,i2]), u([i1,i2]), '-', 'LineWidth', ceil(ss*5), 'Color', col(i,:) );
                  end
               end
            end
            fprintf( '%02d) %d connections\n', p, n );
            export_fig( sprintf( 'fig/symm/%06d.eps', p), '-eps', '-transparent' );
         end
      end

      function shapelets_visualize( this, shapelets )
         NUM_PARTS   = size(shapelets,4);
         NUM_CLASSES = size(shapelets,3);
         for i=1:NUM_PARTS;
            for j=1:NUM_CLASSES;
               I = shapelets(:,:,j,i);
               if sum(I(:)>0)==0; continue; end
               figure;
               imagesc( I );
               title( sprintf( '%s - Part %d', this.clothings(j).name, i ) );
            end
         end
      end

      function [scores] = cpmc_compute( this )
         scores = cell( numel(this.photos_truth), 1 );
         for i=1:numel(this.photos_truth);
            fprintf( '\rcpmc %03d / %03d   ', i, numel(this.photos_truth) );
            pred   = load( sprintf( this.CPMC_PRED,   i ) );
            approx = load( sprintf( this.CPMC_APPROX, i ) );
            weight = pred.potential( pred.objID, 1 );
            mask   = approx.masks( :, :, pred.objID );

            seg    = this.photos_truth(i).segmentation;
            cpmc   = zeros( seg.length, 2 );
            for j=1:seg.length;
               intersection = mask(seg.map==j);
               cpmc(j,1) = weight .* sum(intersection==1) ./ seg.area(j);
               cpmc(j,2) = weight .* sum(intersection==0) ./ seg.area(j);
            end
            scores{i} = cpmc;
         end
         fprintf( '\n' );
      end

      function [X] = unary_features( this )
         X = sbu.ClothParser.feature_transform( this.photos_truth, ...
               'Verbose', true, 'UsePose', true );
      end

      function generate_cpmc_labels( this );
         datapath    = '/share/data/poseseg/';
         outdir      = fullfile(datapath, 'output', 'SegmentationClass29_v2');
         outdir_obj  = fullfile(datapath, 'output', 'SegmentationObject29_v2');

         if ~exist(outdir, 'dir')
             mkdir(outdir);
         end;
         if ~exist(outdir_obj, 'dir')
             mkdir(outdir_obj);
         end;
         %cd(datapath) 
         %startup;

         NUM_CLASSES = numel( this.clothings );
         for i = 1:numel(this.photos_truth);
            fprintf('\r%03d / %03d', i, numel(this.photos_truth));
            ground_truth  = this.photos_truth(i).segmentation.label_map;
            labels        = zeros(size(ground_truth));
            for j=1:NUM_CLASSES;
               labels( ismember( ground_truth, this.clothings(j).id ) ) = j-1;
            end;
            %labels = (ground_truth>0);
            imwrite(uint8(labels), fullfile(outdir,      sprintf('%06d.png', i)))
            imwrite(uint8(labels), fullfile(outdir_obj,  sprintf('%06d.png', i)))
         end;
         fprintf('\n');
      end   

      function generate_cpmc_segments( this )
         %outdir = '/share/data/poseseg/output/cpmccls11/cpmc_segms_150/%06d.mat';
         outdir = 'cpmc_segms_150/%06d.mat';
         for i = 1:numel(this.photos_truth);
            fprintf( '\rcpmc segments %03d / %03d ...   ', i, numel(this.photos_truth) );
            sp       = uint16( this.photos_truth(i).segmentation.map );
            numsp    = max(sp(:));
            sp_app   = logical( diag( ones(numsp,1) ) );
            masks    = zeros( size(sp,1), size(sp,2), numsp );
            for j = 1:numsp;
               masks(:,:,j) = (sp==j);
            end
            masks    = logical( masks );
            save( sprintf(outdir,i), 'sp', 'sp_app', 'masks' );
         end
         fprintf('\n');
      end

      function [scorelist] = compute_L_joint_scores( this )
         joints = this.symm_joints;

         scorelist = cell( numel(this.photos_truth), 1 );
         for p = 1:numel(this.photos_truth);
            fprintf( '\rScores for %03d / %03d...', p, numel(this.photos_truth) );
            seg   = this.photos_truth(p).segmentation;
            uci   = this.photos_truth(p).to_uci;
            box   = pointtobox( uci, this.pa );
            score = zeros( seg.length, size(joints,1), 2 );
            %I     = zeros(size(seg.map));
            for i=1:size(joints,1);
               lmu = uci.point( joints(i,1), : ); 
               rmu = uci.point( joints(i,2), : ); 
               dia = norm( [box.x2(i)-box.x1(i),box.y2(i)-box.y1(i)] );
               sig = this.symm_sigma .* dia .* [1 0;0 1];

               [u,v] = meshgrid( 1:size(seg.map,2), 1:size(seg.map,1) );
               lz    = mvnpdf( [u(:) v(:)], lmu, sig );
               lz    = reshape( lz, size(u,1), size(u,2) );
               rz    = mvnpdf( [u(:) v(:)], rmu, sig );
               rz    = reshape( rz, size(u,1), size(u,2) );
               %I = I + lz + rz;

               for s = 1:seg.length;
                  intersection = (seg.map==s);
                  score(s,i,1) = sum(lz( intersection==1 ));
                  score(s,i,2) = sum(rz( intersection==1 ));
               end
            end
            %figure(1); imshow( I ./ max(I(:)) );
            scorelist{p} = score;
         end
         fprintf( '\n' );
      end

      function [scorelist] = compute_L_limb_scores( this )
         w1    = this.symm_extension;
         w2    = this.symm_width;
         limbs = this.symm_limbs;

         scorelist = cell( numel(this.photos_truth), 1 );
         for p = 1:numel(this.photos_truth);
            fprintf( '\rScores for %03d / %03d...', p, numel(this.photos_truth) );
            seg   = this.photos_truth(p).segmentation;
            uci   = this.photos_truth(p).to_uci;
            score = zeros( seg.length, size(limbs,1), 2 );
            I     = zeros(size(this.photos_truth(p).segmentation.map));
            Il    = zeros(size(this.photos_truth(p).segmentation.map));
            Ir    = zeros(size(this.photos_truth(p).segmentation.map));
            for i=1:size(limbs,1);
               l1 = uci.point( limbs(i,1), : ); 
               l2 = uci.point( limbs(i,2), : ); 
               r1 = uci.point( limbs(i,3), : ); 
               r2 = uci.point( limbs(i,4), : ); 

               ll   = l1-l2;
               lmu  = mean([l1;l2]);
               ltheta = atan2( ll(2), ll(1) );
               c = cos(ltheta); s = sin(ltheta);
               R    = [c -s; s c];
               lsig = R*[w1*(norm(ll)+1) 0; 0 w2]*R';

               rl   = r1-r2;
               rmu  = mean([r1;r2]);
               rtheta = atan2( rl(2), rl(1) );
               c = cos(rtheta); s = sin(rtheta);
               R    = [c -s; s c];
               rsig = R*[w1*(norm(rl)+1) 0; 0 w2]*R';

               [u,v] = meshgrid( 1:size(seg.map,2), 1:size(seg.map,1) );
               lz    = mvnpdf( [u(:) v(:)], lmu, lsig );
               lz    = reshape( lz, size(u,1), size(u,2) );
               rz    = mvnpdf( [u(:) v(:)], rmu, rsig );
               rz    = reshape( rz, size(u,1), size(u,2) );
               Il    = Il + lz;
               Ir    = Ir + rz;
               I     = I + lz + rz;

               for s = 1:seg.length;
                  intersection = (seg.map==s);
                  score(s,i,1) = sum(lz( intersection==1 ));
                  score(s,i,2) = sum(rz( intersection==1 ));
               end
            end
            %clf;
            %II = cat(3, Il, I, Ir );
            %imshow( II ./ max(II(:)) );
            %export_fig( sprintf( 'fig/limbs/%06d.png', p), '-png', '-transparent' );
            scorelist{p} = score;
         end
         fprintf( '\n' );
      end

      function [boxes, photos] = pose_estimation( this )
         photos = this.photos_truth.copy();
         X = photos.to_uci;

         cache = 'fash_final/';
         pe = uci.PoseEstimator( 'name', 'exp_final', ...
               'CachePrefix', cache, 'UseLabels', false );
         pe.train( X(this.trainids), 'FixRNG', 1, ...
               'CachePrefix', cache, 'UseLabels', false );
         pe.model.thresh = min(-2, pe.model.thresh); % Needed for detections
         boxes = pe.estimate( X, ... % Test on _ALL_
               'CachePrefix', cache, 'UseLabels', false );

         % Keep initial pose estimations with associated photos
         for i = 1:numel( photos )
            photos(i).pose = sbu.Pose(boxes{i});
         end
         %save( this.path('pose_estimation.mat'),'S','photos');
      end

      function labels = compute_L_labels( this )
         if this.use_limbs_symm_limbs;
            limbs = this.symm_limbs;
         else
            limbs = this.symm_joints;
         end
         NUM_CLASSES = numel(this.clothings);

         labels = zeros( numel(this.photos_truth), size(limbs,1), 3 );
         for p = 1:numel(this.photos_truth);
            fprintf( '\rlabels %03d / %03d ...', p, numel(this.photos_truth) );
            if this.use_limbs_symm_limbs;
               score = this.X(p).L_score;
            else
               score = this.X(p).J_score;
            end
            seg = this.photos_truth(p).segmentation;
            seglabels = arrayfun( @(s) arrayfun( @(j) ismember( seg.labels(s), this.clothings(j).id ), 1:NUM_CLASSES ), 1:seg.length, 'UniformOutput', false );
            seglabels = cat( 1, seglabels{:} );
            for i = 1:size(limbs,1);
               label1 = sum( seglabels .* repmat( score(:,i,1), 1, size(seglabels,2) ) );
               label2 = sum( seglabels .* repmat( score(:,i,2), 1, size(seglabels,2) ) );
               label3 = label1 + label2;

               [foo,obs1] = max(label1);
               [foo,obs2] = max(label2);
               [foo,obs3] = max(label3);

               labels( p, i, 1 ) = obs1;
               labels( p, i, 2 ) = obs2;
               labels( p, i, 3 ) = obs3;
            end
         end
      end

      function [ logreg, score ] = similarity_visualize( this, photos_similarity )
         if nargin < 2;
            if this.use_real_pose;
               SAVEFILE = this.SIMILARITY_POSE;
            else
               SAVEFILE = this.SIMILARITY;
            end
         end
         data = load( SAVEFILE, 'logreg', 'score' )

         for i=1:numel(data.score);
            seg = this.photos_truth(i).segmentation;
            FNAME = sprintf( 'fig/similarity/%06d.eps', i );
            %if exist( FNAME, 'file' ); continue; end

            score = data.score{i};
            ss = score( (score(:,3) > this.similarity_threshold), : );
            if this.similarity_mintree;
               smin      = ss;
               smin(:,3) = 1-smin(:,3);
               if numel(ss)>0;
                  nMST  = grMinSpanTree( smin );
                  ss    = ss(nMST,:);
               end
            end
            clf;
            seg.show_superpixels();
            hold on;
            for j=1:size(ss,1);
               s1 = ss(j,1);
               s2 = ss(j,2);
               s3 = ss(j,3);

               [I,J] = find( seg.map==s1 );
               u1 = mean(I);
               v1 = mean(J);
               [I,J] = find( seg.map==s2 );
               u2 = mean(I);
               v2 = mean(J);

               plot( [v1,v2], [u1,u2], '-', 'LineWidth', 5, 'Color', 'c' );
            end
            fprintf( '\r%36d) %d connections\n', i, size(ss,1) );
            export_fig( FNAME, '-eps', '-transparent' );
         end
      end

      function [ logreg, score ] = train_similarity( this, photos_similarity )
         if nargin < 2;
            load( this.SIMILARITY_FEATURES, 'photos_similarity' );
            fprintf( 'Loading similarit!\n' );
         end

         trainids = this.trainids;
         testids  = this.testids;
         NPS = numel(this.photos_truth);

         F  = cell( NPS, 1 );
         L  = cell( NPS, 1 );
         score = cell( NPS, 1 );
         for i=1:NPS;
            fprintf( '\rfeatures %03d / %03d', i, NPS );
            p     = this.photos_truth(i);
            seg   = p.segmentation;
            tri   = seg.is_overlapped( p.pose.bbox );
            ps    = photos_similarity{i};
            N     = ((sum(tri)-1)*sum(tri))/2;
            F{i}  = zeros( N, size(ps,3) );
            L{i}  = zeros( N, 1 );
            score{i} = zeros( N, 3 );
            n     = 1;
            for j=1:size(ps,1);
               if ~tri(j); continue; end
               for k=(j+1):size(ps,2);
                  if ~tri(k); continue; end
                  F{i}(n,:) = ps(j,k,:);
                  L{i}(n)   = (seg.labels(j)==seg.labels(k));
                  score{i}(n,1) = j;
                  score{i}(n,2) = k;
                  n = n+1;
               end
            end
         end
         if nargin < 2;
            clear photos_similarity;
         end
         fprintf('\n');
         logreg = liblinear.SVM( cat(1,L{trainids}), cat(1,F{trainids}), ...
               'prob', true, 'type', 0, 'bias', 1, 'autoweight', true );
         for i=1:NPS;
            res            = logreg.predict( cat(1,L{i}), cat(1,F{i}) );
            score{i}(:,3)  = res.val(:,logreg.model.Label==1);
            ids            = (score{i}(:,3) > 0.5); % Filter out bad scores to reduce memory usage
            score{i}       = score{i}(ids,:);
         end

         res = logreg.predict( cat(1,L{testids}), cat(1,F{testids}) );
         out = res.val(:,logreg.model.Label==1);
         Ltrue = cat(1,L{testids});
         ilist = 0.5:0.05:0.9;
         for id=1:numel(ilist);
            i = ilist(id);
            c = confusionmat( logical(cat(1,Ltrue)), (out>i));
            tp = c(1,1);
            fp = c(1,2);
            fn = c(2,1);
            tn = c(2,2);
            precision(id)  = tp / (tp+fp);
            recall(id)     = tp / (tp+fn);
            voc(id)        = tp / (tp+fn+fp);
            connections(id) = sum(out>i);
         end
         ilist
         connections
         precision
         recall
         voc

         if this.use_real_pose;
            SAVEFILE = this.SIMILARITY_POSE;
         else
            SAVEFILE = this.SIMILARITY;
         end
         save( SAVEFILE, 'logreg', 'score' );
      end

      function [ photos_similarity ] = compute_similarity_features( this )
         colourTypes = {'Hsv', 'Lab', 'RGI', 'H', 'Intensity'};

         photos_similarity = cell( numel(this.photos_truth), 1 );
         for i=1:numel(this.photos_truth);
            fprintf( '\rsimilary features %03d / %03d ...', i, numel(this.photos_truth) );

            % Get blobs
            im  = this.photos_truth(i).image;
            seg = this.photos_truth(i).segmentation;
            [blobIndIm blobBoxes neighbours] = seg.blobs;

            % Get image colourspace stuff
            [colourIm imageToSegment]  = Image2ColourSpace(im, colourTypes{1});

            % Get colour histogram
            [colourHist blobSizes]     = BlobStructColourHist (blobIndIm, colourIm); 
            
            % Get texture histogram                                                                         
            textureHist = BlobStructTextureHist(blobIndIm, colourIm);
            % textureHist = BlobStructTextureHistLBP(blobIndIm, colourIm);
           
            % Allocate memory for complete hierarchy.
            numBlobs = size(blobBoxes,1);
            blobStruct.size         = zeros(numBlobs * 2 -1, 1);
            blobStruct.boxes        = zeros(numBlobs * 2 - 1, 4);
            
            % Insert calculated histograms, sizes, and boxes
            blobStruct.size(1:numBlobs)      = blobSizes ./ 3;
            blobStruct.boxes(1:numBlobs,:)   = blobBoxes;
            blobStruct.imSize                = size(im,1) * size(im,2);
            
            % Colours and textures
            blobStruct.textureHist  = zeros(size(textureHist,2), numBlobs * 2 - 1);
            for c = 1:numel(colourTypes);
               col = colourTypes{c};
               % Get image colourspace stuff
               [colourIm imageToSegment]  = Image2ColourSpace(im, col);
               % Get colour histogram
               [colourHist blobSizes]     = BlobStructColourHist(blobIndIm, colourIm); 
               % Get texture histogram                                                                         
               textureHist = BlobStructTextureHist(blobIndIm, colourIm);
               % textureHist = BlobStructTextureHistLBP(blobIndIm, colourIm);
               % Set up histograms
               colHist{c}                 = zeros(size(colourHist, 2), numBlobs * 2 - 1);
               texHist{c}                 = zeros(size(textureHist,2), numBlobs * 2 - 1);
               colHist{c}(:,1:numBlobs)   = colourHist';
               texHist{c}(:,1:numBlobs)   = textureHist';
            end
       
            % Calculate for all
            sim = zeros( seg.length, seg.length, 2+2*numel(colourTypes) );
            for n=1:seg.length;
               for m=(n+1):seg.length;
                  % Size and such
                  box = SSSimBoxFill( n, m, blobStruct );
                  siz = SSSimSize(    n, m, blobStruct );
                  % Put together
                  sim(n,m,1) = box;
                  sim(m,n,1) = box;
                  sim(n,m,2) = siz;
                  sim(m,n,2) = siz;
               end
            end
            for c = 1:numel(colourTypes);
               blobStruct.colourHist   = colHist{c};
               blobStruct.textureHist  = texHist{c};
               for n=1:seg.length;
                  for m=(n+1):seg.length;
                     % Colour and texture
                     col = SSSimColour(  n, m, blobStruct );
                     tex = SSSimTexture( n, m, blobStruct );
                     % Put together
                     id  = 2*c; % two starting are box and sim
                     sim(n,m,id+1) = col;
                     sim(m,n,id+1) = col;
                     sim(n,m,id+2) = tex;
                     sim(m,n,id+2) = tex;
                  end
               end
            end
            photos_similarity{i} = sim;
         end
         fprintf('\n');
         save( this.SIMILARITY_FEATURES, '-v7.3', 'photos_similarity' );
      end


      % Feature saving
      function compute_features( this, path );
         % Compute pose
         fprintf( 'COMPUTING POSE\n' );
         if ~exist([path this.POSE_FILE],'file');
            boxes = this.pose_estimation();
            save( [path this.POSE_FILE], 'boxes' );
         end

         % Get features
         fprintf( 'COMPUTING FEATURES\n' );
         if ~exist([path this.FEATURES_FILE],'file');
            this.use_real_pose = 0;
            this = this.features_load;
            Xtruth = this.unary_features();
            save( [path this.FEATURES_FILE], 'Xtruth' );
         end
         if ~exist([path this.FEATURES_POSE_FILE],'file');
            this.use_real_pose = 1;
            this = this.features_load;
            Xtruth = this.unary_features();
            save( [path this.FEATURES_POSE_FILE], 'Xtruth' );
         end
         this = this.features_load; % Reload the features

         % Tamara's logreg
         fprintf( 'COMPUTING TAMARA LOGREGS\n' );
         this.use_real_pose = 0;
         this.compute_logreg_unaries( path );
         this.use_real_pose = 1;
         this.compute_logreg_unaries( path );

         % Similarity
         [photos_similarity] = this.compute_similarity_features();
         this.use_real_pose = 0;
         this = this.features_load;
         this.train_similarity( photos_similarity );
         this.use_real_pose = 1;
         this = this.features_load;
         this.train_similarity( photos_similarity );
         clear photos_similarity;

         % Shapelets
         fprintf( 'COMPUTING SHAPELETS\n' );
         if ~exist([path this.SHAPELETS_FILE],'file') || ~exist([path this.SHAPELETS_POSE_FILE]);
            this.use_real_pose = 0;
            this        = this.features_load();
            shapelets   = this.shapelets_compute();
            if ~exist([path this.SHAPELETS_FILE],'file');
               score       = this.shapelets_evaluate( this.photos_truth, shapelets );
               save( [path this.SHAPELETS_FILE], 'score' );
            end
            if ~exist([path this.SHAPELETS_POSE_FILE],'file');
               this.use_real_pose = 1;
               this        = this.features_load();
               score       = this.shapelets_evaluate( this.photos_truth, shapelets );
               save( [path this.SHAPELETS_POSE_FILE], 'score' );
            end
         end

         % CPMC Unaries - ignores pose
         fprintf( 'COMPUTING CPMC LOGREGS\n' );
         this.compute_cpmc_logreg_unaries( path );

         % Independent of pose
         fprintf( 'COMPUTING CPMC\n' );
         if ~exist([path this.CPMC_FILE],'file');
            score = this.cpmc_compute();
            save( [path this.CPMC_FILE], 'score' );
         end
         %score = e.cpmc_compute(); save( [e.FEATPATH e.CPMC_FILE], 'score' );
      end
   end
end




