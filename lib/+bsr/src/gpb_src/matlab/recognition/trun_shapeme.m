%K = 1000;
K = 2857;
[train, test] = create_train_test('/home/mmaire_data/caltech101_gb_pb_pca','/home/mmaire_data/caltech101_resized',15,15);
train_shape = train;

%train_shape.filenames     = train_shape.filenames(1:15:end);
%train_shape.filenames_img = train_shape.filenames_img(1:15:end);
%train_shape.class_id      = train_shape.class_id(1:15:end);
%train_shape.model_id      = train_shape.model_id(1:15:end);
%train_shape.n_train_per_class = 1;

train_shape.filenames     = train_shape.filenames(1:150:end);
train_shape.filenames_img = train_shape.filenames_img(1:150:end);
train_shape.class_id      = train_shape.class_id(1:150:end);
train_shape.model_id      = train_shape.model_id(1:150:end);
train_shape.n_train_per_class = 1;

train = train_shape;
test = train_shape;
test.n_test_per_class = 11/101;

build_vocab_shapeme(train_shape,'pca');
build_db_shapeme(train,'pca',K);
[cmx, score] = run_test_shapeme(train,test,'pca');
