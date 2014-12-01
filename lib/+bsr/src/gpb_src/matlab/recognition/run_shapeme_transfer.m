K = 1000;
[train, test] = create_train_test('/home/mmaire_data/caltech101_gb_pb_pca','/home/mmaire_data/caltech101_resized',15,15);

% get training set for shapemes
shape_inds = 751:15:1515;
train_shape = train;
train_shape.filenames     = train_shape.filenames(shape_inds);
train_shape.filenames_img = train_shape.filenames_img(shape_inds);
train_shape.class_id      = train_shape.class_id(shape_inds);
train_shape.model_id      = train_shape.model_id(shape_inds);
train_shape.n_train_per_class = 1;

% get training set for recognition
rec_inds = 1:750;
train.filenames     = train.filenames(rec_inds);
train.filenames_img = train.filenames_img(rec_inds);
train.class_id      = train.class_id(rec_inds);
train.model_id      = train.model_id(rec_inds);

% get test set for recognition
test.filenames     = test.filenames(rec_inds);
test.filenames_img = test.filenames_img(rec_inds);
test.class_id      = test.class_id(rec_inds);
test.model_id      = test.model_id(rec_inds);

% run experiment
build_vocab_shapeme(train_shape,'pca');
build_db_shapeme(train,'pca',K);
[cmx, score] = run_test_shapeme(train,test,'pca');
